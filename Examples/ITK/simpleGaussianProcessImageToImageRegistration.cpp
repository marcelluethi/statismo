/*
 * This file is part of the statismo library.
 *
 * Author: Christoph Jud (christoph.jud@unibas.ch)
 *
 * Copyright (c) 2011 University of Basel
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the project's author nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * 
 * This simple example is to illustrate how a Gaussian Process model can be used for image registration.
 * For a more complete example, see the itkLowRankGaussianProcessImageToImageRegistration example in the same folder.
 *
 * Usage:
 * ./bin/itkSimpleGaussianProcessImageToImageRegistration share/data/hand_images/hand-1.vtk share/data/hand_images/hand-2.vtk /tmp/deformationfield.vtk
 *
 */

#include "itkVectorImageRepresenter.h"
#include "statismo_ITK/itkStatisticalModel.h"
#include "statismo_ITK/itkLowRankGPModelBuilder.h"
#include "statismo_ITK/itkDataManager.h"
#include "statismo_ITK/itkInterpolatingStatisticalDeformationModelTransform.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkLBFGSOptimizer.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageRegistrationMethod.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCommand.h"
#include "itkImage.h"

#include "itkDirectory.h"
#include <sys/types.h>
#include <errno.h>
#include <iostream>
#include <iomanip>
#include <string>

#define GaussianSigma 70
#define GaussianScale 100
#define NumBasisFunctions 100
#define NumIterations 100

const unsigned int Dimension = 2;
typedef itk::Image<float, Dimension> ImageType;
typedef itk::Image<itk::Vector<float, Dimension>, Dimension> VectorImageType;

typedef itk::VectorImageRepresenter<float, Dimension, Dimension> RepresenterType;
typedef itk::StatisticalModel<RepresenterType> StatisticalModelType;
typedef itk::InterpolatingStatisticalDeformationModelTransform<RepresenterType,
		double, Dimension> TransformType;

typedef itk::ImageFileReader<ImageType> ImageReaderType;
typedef itk::MeanSquaresImageToImageMetric<ImageType, ImageType> MetricType;
typedef itk::LBFGSOptimizer OptimizerType;
typedef itk::ImageRegistrationMethod<ImageType, ImageType> RegistrationFilterType;
typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;



/**
 * A scalar valued gaussian kernel.
 */
template <class TPoint>
class GaussianKernel: public statismo::ScalarValuedKernel<TPoint>{
public:
	typedef typename  TPoint::CoordRepType CoordRepType;
	typedef vnl_vector<CoordRepType> VectorType;

	GaussianKernel(double sigma) : m_sigma(sigma), m_sigma2(sigma * sigma) {}

	inline double operator()(const TPoint& x, const TPoint& y) const {
		VectorType xv = x.GetVnlVector();
		VectorType yv = y.GetVnlVector();

		VectorType r = yv - xv;
		return exp(-dot_product(r, r) / m_sigma2);
	}

	std::string GetKernelInfo() const {
		std::ostringstream os;
		os << "GaussianKernel(" << m_sigma << ")";
		return os.str();
	}

private:

	double m_sigma;
	double m_sigma2;
};




/*
 * Build a low-rank Gaussian process model using a Gaussian kernel function
 * 
 * Input:
 *  - Filename to reference image
 *  
 * Output:
 *  - a (statismo) statistical model.
 *  
 * Comment:
 *  - The composition of different kernel functions can be even more complicated. For example
 *    a linear combination of different kernel functions is again a kernel function and thus
 *    can be handled by the LowRankGPModelBuilder.
 */
StatisticalModelType::Pointer buildLowRankGPModel(
		const char* referenceFilename) {

	typedef itk::LowRankGPModelBuilder<RepresenterType> ModelBuilderType;
	typedef std::vector<std::string> StringVectorType;
	typedef itk::ImageFileReader<VectorImageType> DeformationFieldReaderType;
	typedef RepresenterType::PointType PointType;

	std::cout << "Building low-rank Gaussian process deformation model... "
			<< std::flush;

	// we take an arbitrary dataset as the reference, as they have all the same resolution anyway
	DeformationFieldReaderType::Pointer referenceReader =
			DeformationFieldReaderType::New();
	referenceReader->SetFileName(referenceFilename);
	referenceReader->Update();

	RepresenterType::Pointer representer = RepresenterType::New();
	representer->SetReference(referenceReader->GetOutput());

	const GaussianKernel<PointType> gk = GaussianKernel<PointType>(GaussianSigma); // a Gaussian kernel with sigma=gaussianKernelSigma

	// make the kernel matrix valued and scale it by a factor of 100
	const statismo::MatrixValuedKernel<PointType>& mvGk = statismo::UncorrelatedMatrixValuedKernel<PointType>(&gk, representer->GetDimensions());
	const statismo::MatrixValuedKernel<PointType>& scaledGk = statismo::ScaledKernel<PointType>(&mvGk, GaussianScale); // apply Gaussian scale parameter

	ModelBuilderType::Pointer gpModelBuilder = ModelBuilderType::New();
	gpModelBuilder->SetRepresenter(representer);
	StatisticalModelType::Pointer model = gpModelBuilder->BuildNewZeroMeanModel(scaledGk, NumBasisFunctions); // number of basis functions

	std::cout << "[done]" << std::endl;
	return model;
}

/*
 * Image to image registration method using a statismo statistical model.
 * 
 * The standard parametric registration framework of ITK is used for registration, where
 * the transform is a InterpolatingStatisticalDeformationModelTransform.
 * 
 * Input:
 *  - Filename of reference image.
 *  - Filename of target image.
 *  - Statistical model.
 *  
 * Output:
 *  - The deformation field.
 */
VectorImageType::Pointer modelBasedImageToImageRegistration(
		std::string referenceFilename, std::string targetFilename,
		StatisticalModelType::Pointer model) {

	ImageReaderType::Pointer referenceReader = ImageReaderType::New();
	referenceReader->SetFileName(referenceFilename.c_str());
	referenceReader->Update();
	ImageType::Pointer referenceImage = referenceReader->GetOutput();

	ImageReaderType::Pointer targetReader = ImageReaderType::New();
	targetReader->SetFileName(targetFilename.c_str());
	targetReader->Update();
	ImageType::Pointer targetImage = targetReader->GetOutput();

	// do the fitting
	TransformType::Pointer transform = TransformType::New();
	transform->SetStatisticalModel(model);
	transform->SetIdentity();

	// Setting up the fitting
	OptimizerType::Pointer optimizer = OptimizerType::New();
	optimizer->MinimizeOn();
	optimizer->SetMaximumNumberOfFunctionEvaluations(NumIterations);

	MetricType::Pointer metric = MetricType::New();
	InterpolatorType::Pointer interpolator = InterpolatorType::New();

	RegistrationFilterType::Pointer registration =
			RegistrationFilterType::New();
	registration->SetInitialTransformParameters(transform->GetParameters());
	registration->SetMetric(metric);
	registration->SetOptimizer(optimizer);
	registration->SetTransform(transform);
	registration->SetInterpolator(interpolator);
	registration->SetFixedImage(targetImage);
	registration->SetFixedImageRegion(targetImage->GetBufferedRegion());
	registration->SetMovingImage(referenceImage);

	try {
		std::cout << "Performing registration... " << std::flush;
		registration->Update();
		std::cout << "[done]" << std::endl;

	} catch (itk::ExceptionObject& o) {
		std::cout << "caught exception " << o << std::endl;
	}

	return model->DrawSample(transform->GetCoefficients());
}

int main(int argc, char* argv[]) {

	std::string referenceFilename("");
	std::string targetFilename("");
	std::string outputFilename("");

	// parse command line parameters
	if (argc != 4) {
		std::cout << "usage\t" << argv[0]
				<< " referenceFilename targetFilename outputFilename"
				<< std::endl;
		exit(-1);
	}

	referenceFilename = argv[1];
	targetFilename = argv[2];
	outputFilename = argv[3];

	// perform low-rank approximation of Gaussian process prior
	StatisticalModelType::Pointer model = buildLowRankGPModel(
			referenceFilename.c_str());

	// perform image to image registration using the Gaussian process deformation model
	VectorImageType::Pointer deformationField =
			modelBasedImageToImageRegistration(referenceFilename,
					targetFilename, model);

	// write deformation field
	itk::ImageFileWriter<VectorImageType>::Pointer df_writer =
			itk::ImageFileWriter < VectorImageType > ::New();
	df_writer->SetFileName(outputFilename);
	df_writer->SetInput(deformationField);
	df_writer->Update();

	std::cout
			<< "Low-rank Gaussian process image to image registration has been successfully finished."
			<< std::endl;
}

