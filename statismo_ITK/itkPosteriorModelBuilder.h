/*
 * This file is part of the statismo library.
 *
 * Author: Marcel Luethi (marcel.luethi@unibas.ch)
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
 */


#ifndef ITK_POSTERIOR_MODELBUILDER_H_
#define ITK_POSTERIOR_MODELBUILDER_H_

#include "itkObject.h"
#include "statismoITKConfig.h"
#include "itkDataManager.h"
#include "itkStatisticalModel.h"
#include "statismo/PosteriorModelBuilder.h"


namespace itk
{

/**
 * \brief ITK Wrapper for the statismo::PosteriorModelBuilder class.
 * \see statismo::PosteriorModelBuilder for detailed documentation.
 */
template <class Representer>
class PosteriorModelBuilder : public Object {
public:

	typedef PosteriorModelBuilder            Self;
	typedef Object	Superclass;
	typedef SmartPointer<Self>                Pointer;
	typedef SmartPointer<const Self>          ConstPointer;

	itkNewMacro( Self );
	itkTypeMacro( PosteriorModelBuilder, Object );

	typedef statismo::PosteriorModelBuilder<Representer> ImplType;
	typedef statismo::DataManager<Representer> DataManagerType;
	typedef typename DataManagerType::SampleDataStructureListType SampleDataStructureListType;



	template <class F>
	typename std::tr1::result_of<F()>::type callstatismoImpl(F f) const {
		try {
			  return f();
		}
		 catch (statismo::StatisticalModelException& s) {
			itkExceptionMacro(<< s.what());
		}
	}


	PosteriorModelBuilder() : m_impl(ImplType::Create()) {}

	virtual ~PosteriorModelBuilder() {
		if (m_impl) {
			delete m_impl;
			m_impl = 0;
		}
	}


	// create statismo stuff
	typedef typename Representer::ValueType ValueType;
	typedef typename Representer::PointType PointType;
	typedef typename statismo::PosteriorModelBuilder<Representer>::PointValueListType PointValueListType;
	typedef typename statismo::PosteriorModelBuilder<Representer>::PointValueWithCovariancePairType PointValueWithCovariancePairType;
	typedef typename statismo::PosteriorModelBuilder<Representer>::PointValueWithCovarianceListType PointValueWithCovarianceListType;



	typename StatisticalModel<Representer>::Pointer BuildNewModelFromModel(const StatisticalModel<Representer>* model, const PointValueListType& pointValues, double pointValuesNoiseVariance,  bool computeScores=true) {
		statismo::StatisticalModel<Representer>* model_statismo = model->GetstatismoImplObj();
		statismo::StatisticalModel<Representer>* new_model_statismo = callstatismoImpl(std::tr1::bind(&ImplType::BuildNewModelFromModel, this->m_impl, model_statismo, pointValues, pointValuesNoiseVariance, computeScores));
		typename StatisticalModel<Representer>::Pointer model_itk = StatisticalModel<Representer>::New();
		model_itk->SetstatismoImplObj(new_model_statismo);
		return model_itk;
	}

	typename StatisticalModel<Representer>::Pointer BuildNewModel(SampleDataStructureListType SampleDataStructureList, const PointValueListType& pointValues, double pointValuesNoiseVariance, double noiseVariance) {
		statismo::StatisticalModel<Representer>* model_statismo = callstatismoImpl(std::tr1::bind(&ImplType::BuildNewModel, this->m_impl, SampleDataStructureList ,pointValues, pointValuesNoiseVariance, noiseVariance));
		typename StatisticalModel<Representer>::Pointer model_itk = StatisticalModel<Representer>::New();
		model_itk->SetstatismoImplObj(model_statismo);
		return model_itk;
	}

	typename StatisticalModel<Representer>::Pointer BuildNewModelFromModel(const StatisticalModel<Representer>* model, const PointValueWithCovarianceListType& pointValuesWithCovariance,  bool computeScores=true) {
		statismo::StatisticalModel<Representer>* model_statismo = model->GetstatismoImplObj();
		statismo::StatisticalModel<Representer>* new_model_statismo = callstatismoImpl(std::tr1::bind(&ImplType::BuildNewModelFromModel, this->m_impl, model_statismo, pointValuesWithCovariance, computeScores));
		typename StatisticalModel<Representer>::Pointer model_itk = StatisticalModel<Representer>::New();
		model_itk->SetstatismoImplObj(new_model_statismo);
		return model_itk;
	}

	typename StatisticalModel<Representer>::Pointer BuildNewModel(const SampleDataStructureListType& SampleDataStructureList, const PointValueWithCovarianceListType& pointValuesWithCovariance, 	double noiseVariance) {
		statismo::StatisticalModel<Representer>* model_statismo = callstatismoImpl(std::tr1::bind(&ImplType::BuildNewModel, this->m_impl, SampleDataStructureList, pointValuesWithCovariance, noiseVariance));
		typename StatisticalModel<Representer>::Pointer model_itk = StatisticalModel<Representer>::New();
		model_itk->SetstatismoImplObj(model_statismo);
		return model_itk;
	}

private:
	PosteriorModelBuilder(const PosteriorModelBuilder& orig);
	PosteriorModelBuilder& operator=(const PosteriorModelBuilder& rhs);

	ImplType* m_impl;
};


}

#endif /* ITK_POSTERIOR_MODEL_BUILDER */
