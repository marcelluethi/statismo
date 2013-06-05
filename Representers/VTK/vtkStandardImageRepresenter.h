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



#ifndef VTKSTANDARDIMAGERESENTER_H_
#define VTKSTANDARDIMAGERESENTER_H_


#include "vtkStructuredPoints.h"
#include "vtkPoint.h"
#include "vtkPixel.h"
#include "statismo/CommonTypes.h"
#include "statismo/Domain.h"
#include "vtkSmartPointer.h"

#include <H5Cpp.h>

/**
 * \brief Representer class for vtkStructuredPoints of arbitrary scalar type and PixelDimension
 * \sa Representer
 */



template <class TScalar, unsigned PixelDimensions>
class vtkStandardImageRepresenter  {
public:



	typedef vtkStructuredPoints DatasetType;
	typedef vtkStructuredPoints* DatasetPointerType;
	typedef const vtkStructuredPoints* DatasetConstPointerType;

	//typedef double* PointType;
	typedef vtkPoint PointType;

	typedef statismo::Domain<PointType> DomainType;

	typedef vtkNDPixel<TScalar, PixelDimensions> ValueType;

	static vtkStandardImageRepresenter* Create(DatasetConstPointerType reference) { return new vtkStandardImageRepresenter(reference); }

	static vtkStandardImageRepresenter* Load(const H5::CommonFG& fg);
	vtkStandardImageRepresenter* Clone() const;

	struct DatasetInfo {}; // not used for this representer, but needs to be here as it is part of the generic interface

	virtual ~vtkStandardImageRepresenter();
	void Delete() const { delete this; }



	static unsigned GetDimensions() { return  PixelDimensions; }
	const DomainType& GetDomain() const  { return m_domain; }

	static std::string GetName() { return "vtkStandardImageRepresenter"; }

	const vtkStructuredPoints* GetReference() const { return m_reference; }


	DatasetPointerType DatasetToSample(DatasetConstPointerType ds, DatasetInfo* notUsed) const;
	statismo::VectorType SampleToSampleVector(DatasetConstPointerType sample) const;
	DatasetPointerType SampleVectorToSample(const statismo::VectorType& sample) const;


	ValueType PointSampleFromSample(DatasetConstPointerType sample, unsigned ptid) const;
	statismo::VectorType PointSampleToPointSampleVector(const ValueType& v) const;
	ValueType PointSampleVectorToPointSample(const statismo::VectorType& samplePoint) const;

	unsigned GetPointIdForPoint(const PointType& pt) const;

	unsigned GetNumberOfPoints() const;
	void Save(const H5::CommonFG& fg) const;

	static void DeleteDataset(DatasetPointerType  d);
    static unsigned GetNumberOfPoints(DatasetPointerType  reference);

	 /* Maps a (Pointid,component) tuple to a component of the internal matrix.
	 * This is used to locate the place in the matrix to store the elements for a given point.
	 * @params ptId The point id
	 * @params the Component Index (range 0, Dimensionality)
	 * @returns an index.
	 */
	static unsigned MapPointIdToInternalIdx(unsigned ptId, unsigned componentInd) {
		return ptId * GetDimensions() + componentInd;
	}

private:

	vtkStandardImageRepresenter(DatasetConstPointerType reference);

	DatasetType* m_reference;
	DomainType m_domain;
};

#include "vtkStandardImageRepresenter.txx"

#endif /* VTKStandardImageRepRESENTER_H_ */
