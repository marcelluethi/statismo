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
 
%module statismo

%include "typemaps.i"
%include "std_string.i"
%include "std_list.i"
%include "std_vector.i"
%include "std_pair.i"
%include "std_vector.i"
%include "carrays.i"

%include "statismoTypemaps.i"
%include "vtkPolyDataRepresenter.i"
%include "vtkUnstructuredGridRepresenter.i"
%include "vtkStructuredPointsRepresenter.i"
%include "TrivialVectorialRepresenter.i"

  
%{
#include "statismo/DataManager.h"
#include "statismo/DataManagerWithSurrogates.h"
#include "statismo/StatisticalModel.h"
#include "statismo/PartiallyFixedModelBuilder.h"
#include "statismo/PosteriorModelBuilder.h"
#include "statismo/ReducedVarianceModelBuilder.h"
#include "statismo/PCAModelBuilder.h"
#include "statismo/Exceptions.h"
#include "statismo/CommonTypes.h"
#include <list>
#include <string>
%}


namespace statismo { 
class StatisticalModelException {
public:
	StatisticalModelException(const char* message);
};
}

%exception {
   try {
      $action
   } catch (statismo::StatisticalModelException &e) {
      PyErr_SetString(PyExc_RuntimeError, const_cast<char*>(e.what()));
      return NULL;
   }
}



//////////////////////////////////////////////////////
// DataManager
//////////////////////////////////////////////////////

namespace statismo { 
template <typename Representer>
struct SampleDataStructure {
	typedef typename Representer::DatasetPointerType DatasetPointerType;
	
	std::string GetDatasetURI() const;
	%newobject GetSample;
	const DatasetPointerType GetSample() const;

	const  statismo::VectorType& GetSampleVector() const;
	virtual ~SampleDataStructure();
private:
	SampleDataStructure(const Representer* representer, const std::string& filename, const VectorType& sampleVector);
};

template <typename Representer>
struct SampleDataStructureWithSurrogates : public SampleDataStructure<Representer> { 
	const VectorType& GetSurrogateVector() const;
	const std::string& GetSurrogateFilename() const;
};
}

%template(SampleDataStructure_tvr) statismo::SampleDataStructure<TrivialVectorialRepresenter>;
%template(SampleDataStructure_vtkPD) statismo::SampleDataStructure<vtkPolyDataRepresenter>;
%template(SampleDataStructure_vtkUG) statismo::SampleDataStructure<vtkUnstructuredGridRepresenter>;
%template(SampleDataStructure_vtkSPF3) statismo::SampleDataStructure<vtkStructuredPointsRepresenter<float, 3> >;
%template(SampleDataStructure_vtkSPSS1) statismo::SampleDataStructure<vtkStructuredPointsRepresenter<signed short, 1> >;
%template(SampleDataStructureWithSurrogates_tvr) statismo::SampleDataStructureWithSurrogates<TrivialVectorialRepresenter>;
%template(SampleDataStructureWithSurrogates_vtkPD) statismo::SampleDataStructureWithSurrogates<vtkPolyDataRepresenter>;
%template(SampleDataStructureWithSurrogates_vtkUG) statismo::SampleDataStructureWithSurrogates<vtkUnstructuredGridRepresenter>;
%template(SampleDataStructureWithSurrogates_vtkSPF3) statismo::SampleDataStructureWithSurrogates<vtkStructuredPointsRepresenter<float, 3> >;
%template(SampleDataStructureWithSurrogates_vtkSPSS1) statismo::SampleDataStructureWithSurrogates<vtkStructuredPointsRepresenter<signed short, 1> >;


%traits_swigtype(statismo::SampleDataStructure<vtkPolyDataRepresenter>); // workaround for a swig bug with const ptr
%fragment(SWIG_Traits_frag(statismo::SampleDataStructure<vtkPolyDataRepresenter>));  // workaround for a swig bug with const ptr
%template(SampleDataStructureList_vtkPD) std::list<const statismo::SampleDataStructure<vtkPolyDataRepresenter> *>;

%traits_swigtype(statismo::SampleDataStructure<vtkUnstructuredGridRepresenter>); // workaround for a swig bug with const ptr
%fragment(SWIG_Traits_frag(statismo::SampleDataStructure<vtkUnstructuredGridRepresenter>));  // workaround for a swig bug with const ptr
%template(SampleDataStructureList_vtkUG) std::list<const statismo::SampleDataStructure<vtkUnstructuredGridRepresenter> *>;

%traits_swigtype(statismo::SampleDataStructure<vtkStructuredPointsRepresenter<float, 3> >);  // workaround for a swig bug with const ptr
%fragment(SWIG_Traits_frag(statismo::SampleDataStructure<vtkStructuredPointsRepresenter<float, 3> >));  // workaround for a swig bug with const ptr
%template(SampleDataStructureList_vtkSPF3) std::list<const statismo::SampleDataStructure<vtkStructuredPointsRepresenter<float, 3> >* > ;

%traits_swigtype(statismo::SampleDataStructure<vtkStructuredPointsRepresenter<signed short, 1> >);  // workaround for a swig bug with const ptr
%fragment(SWIG_Traits_frag(statismo::SampleDataStructure<vtkStructuredPointsRepresenter<signed short, 1> >));  // workaround for a swig bug with const ptr
%template(SampleDataStructureList_vtkSPSS1) std::list<const statismo::SampleDataStructure<vtkStructuredPointsRepresenter<signed short, 1> >* > ;


%traits_swigtype(statismo::SampleDataStructure<TrivialVectorialRepresenter>);  // workaround for a swig bug with const ptr
%fragment(SWIG_Traits_frag(statismo::SampleDataStructure<TrivialVectorialRepresenter>));  // workaround for a swig bug with const ptr
%template(SampleDataStructureList_tvr) std::list<const statismo::SampleDataStructure<TrivialVectorialRepresenter>* > ;


template <typename Representer>
class CrossValidationFold {
public:
	typedef SampleDataStructure<Representer> SampleDataStructureType;
	typedef std::list<const SampleDataStructureType*> SampleDataStructureListType;

	CrossValidationFold();
	SampleDataStructureListType GetTrainingData() const;
	SampleDataStructureListType GetTestingData() const;

};
%template(CrossValidationFold_tvr) statismo::CrossValidationFold<TrivialVectorialRepresenter>;
%template(CrossValidationFold_vtkPD) statismo::CrossValidationFold<vtkPolyDataRepresenter>;
%template(CrossValidationFold_vtkUG) statismo::CrossValidationFold<vtkUnstructuredGridRepresenter>;
%template(CrossValidationFold_vtkSPF3) statismo::CrossValidationFold<vtkStructuredPointsRepresenter<float, 3> >;
%template(CrossValidationFold_vtkSPSS1) statismo::CrossValidationFold<vtkStructuredPointsRepresenter<signed short, 1> >;

%template(CrossValidationFoldList_tvr) std::list<statismo::CrossValidationFold<TrivialVectorialRepresenter> >;
%template(CrossValidationFoldList_vtkPD) std::list<statismo::CrossValidationFold<vtkPolyDataRepresenter> >;
%template(CrossValidationFoldList_vtkUG) std::list<statismo::CrossValidationFold<vtkUnstructuredGridRepresenter> >;
%template(CrossValidationFoldList_vtkSPF3) std::list<statismo::CrossValidationFold<vtkStructuredPointsRepresenter<float, 3> > > ;
%template(CrossValidationFoldList_vtkSPSS1) std::list<statismo::CrossValidationFold<vtkStructuredPointsRepresenter<signed short, 1> > > ;


//////////////////////////////////////////////////////
// DataManager
//////////////////////////////////////////////////////

namespace statismo {
template <typename Representer>
class DataManager {
public:

	typedef SampleDataStructure<Representer> SampleDataStructureType;
	typedef std::list<const SampleDataStructureType *> SampleDataStructureListType;

	typedef CrossValidationFold<Representer> CrossValidationFoldType;
	typedef std::list<CrossValidationFoldType> CrossValidationFoldListType;



	
	%newobject Create;
	static DataManager* Create(const Representer*);
	static DataManager* Load(const char* filename);
		
	void AddDataset(Representer::DatasetConstPointerType dataset, const std::string& URI);
	unsigned GetNumberOfSamples() const;

	void Save(const char* filename);

	
	/** cross validation functionality */
	CrossValidationFoldListType GetCrossValidationFolds(unsigned nFolds, bool randomize = true) const;
	SampleDataStructureListType GetSampleDataStructure() const;
	
private:
	DataManager();
};
}
%template(DataManager_tvr) statismo::DataManager<TrivialVectorialRepresenter>;
%template(DataManager_vtkPD) statismo::DataManager<vtkPolyDataRepresenter>;
%template(DataManager_vtkUG) statismo::DataManager<vtkUnstructuredGridRepresenter>;
%template(DataManager_vtkSPF3) statismo::DataManager<vtkStructuredPointsRepresenter<float, 3> >;
%template(DataManager_vtkSPSS1) statismo::DataManager<vtkStructuredPointsRepresenter<signed short, 1> >;


namespace statismo {
template <typename Representer>
class DataManagerWithSurrogates : public DataManager<Representer> {
public:


	static DataManagerWithSurrogates<Representer>* Create(const Representer* representer, const std::string& surrogTypeFilename);
	void AddDatasetWithSurrogates(Representer::DatasetConstPointerType datasetFilename, const std::string& surrogateFilename,  const std::string& surrogateFilename);
		
private:
	DataManagerWithSurrogates();
};
}
%template(DataManagerWithSurrogates_tvr) statismo::DataManagerWithSurrogates<TrivialVectorialRepresenter>;
%template(DataManagerWithSurrogates_vtkPD) statismo::DataManagerWithSurrogates<vtkPolyDataRepresenter>;
%template(DataManagerWithSurrogates_vtkUG) statismo::DataManagerWithSurrogates<vtkUnstructuredGridRepresenter>;
%template(DataManagerWithSurrogates_vtkSPF3) statismo::DataManagerWithSurrogates<vtkStructuredPointsRepresenter<float, 3> >;
%template(DataManagerWithSurrogates_vtkSPSS1) statismo::DataManagerWithSurrogates<vtkStructuredPointsRepresenter<signed short, 1> >;


//////////////////////////////////////////////////////
// PointValuePair
//////////////////////////////////////////////////////


%template(PointValuePair_tvr) std::pair<TrivialVectorialRepresenter::PointType, TrivialVectorialRepresenter::ValueType>;
%template(PointValuePair_vtkPD) std::pair<vtkPolyDataRepresenter::PointType, vtkPolyDataRepresenter::ValueType>;
%template(PointValueList_vtkPD) std::list<std::pair<vtkPolyDataRepresenter::PointType, vtkPolyDataRepresenter::ValueType> >;
//%template(PointValuePair_vtkUG) std::pair<vtkUnstructuredGridRepresenter::PointType, vtkUnstructuredGridRepresenter::ValueType>;
//%template(PointValueList_vtkUG) std::list<std::pair<vtkUnstructuredGridRepresenter::PointType, vtkUnstructuredGridRepresenter::ValueType> >;
%template(PointIdValuePair_tvr) std::pair<unsigned, TrivialVectorialRepresenter::ValueType>;
%template(PointIdValuePair_vtkPD) std::pair<unsigned, vtkPolyDataRepresenter::ValueType>;
%template(PointIdValueList_vtkPD) std::list<std::pair<unsigned, vtkPolyDataRepresenter::ValueType> >;
//%template(PointIdValuePair_vtkUG) std::pair<unsigned, vtkUnstructuredGridRepresenter::ValueType>;
//%template(PointIdValueList_vtkUG) std::list<std::pair<unsigned, vtkUnstructuredGridRepresenter::ValueType> >;

//////////////////////////////////////////////////////
// PointValueWithCovariancePair
//////////////////////////////////////////////////////

// Never got this to work properly. The wrapping of statismo::MatrixType is not sophisticated enough to work with the std::pair and list wrappers.

//%template(PointValuePair_tvr) std::pair<TrivialVectorialRepresenter::PointType, TrivialVectorialRepresenter::ValueType>;
//%template(PointValueWithCovariancePair_vtkPD) std::pair< std::pair< vtkPoint,vtkPoint >,statismo::MatrixType >; 
//%template(PointValueWithCovarianceList_vtkPD) std::list< std::pair< std::pair< vtkPoint,vtkPoint >,statismo::MatrixType > >;
//%template(MatrixList) std::list<statismo::MatrixType>;
//%template(MatrixMatrixPair) std::pair< statismo::MatrixType, statismo::MatrixType >;
//%template(MatrixMatrixList) std::list< std::pair< statismo::MatrixType, statismo::MatrixType > >;
//%template(PointMatrixPair) std::pair< vtkPoint, statismo::MatrixType >;
//%template(PointMatrixList) std::list< std::pair< vtkPoint, statismo::MatrixType > >;


//////////////////////////////////////////////////////
// ModelInfo
//////////////////////////////////////////////////////

namespace statismo { 

class BuilderInfo {
public:
	typedef std::pair<std::string, std::string> KeyValuePair;
	typedef std::list<KeyValuePair> KeyValueList;
			
	const KeyValueList& GetDataInfo() const;
	const KeyValueList& GetParameterInfo() const;	
};

class ModelInfo {
public:
	typedef std::vector<BuilderInfo> BuilderInfoList;
	const statismo::MatrixType& GetScoresMatrix() const;

	virtual void Save(const H5::CommonFG& publicFg) const;
	virtual void Load(const H5::CommonFG& publicFg);			
	BuilderInfoList GetBuilderInfoList() const ;
};
}
%template (KeyValuePair) std::pair<std::string, std::string>; 
%template(KeyValueList) std::list<std::pair<std::string, std::string> >;
%template(BuilderInfoList) std::vector<statismo::BuilderInfo>;

/////////////////////////////////////////////////////////////////
// Domain
/////////////////////////////////////////////////////////////////
namespace statismo { 
template <typename PointType>
class Domain {
public:
	typedef std::vector<PointType> DomainPointsListType;
	const DomainPointsListType& GetDomainPoints() const;
	const unsigned GetNumberOfPoints() const;
};
}
%template(DomainVtkPoint) statismo::Domain<vtkPoint>;
%template(DomainId) statismo::Domain<unsigned int>;
%template(DomainPointsListVtkPoint) std::vector<vtkPoint>;
%template(DomainPointsListId) std::vector<unsigned int>;

//////////////////////////////////////////////////////
// StatisticalModel
//////////////////////////////////////////////////////

namespace statismo { 
template <class Representer>
class StatisticalModel {
public:

     typedef Representer::DatasetPointerType DatasetPointerType;
     typedef Representer::DatasetConstPointerType DatasetConstPointerType;     

	typedef  std::pair<typename Representer::PointType, typename Representer::ValueType>  PointValuePairType;
 	typedef std::list<PointValuePairType> PointValueListType;
	typedef  std::pair<unsigned, typename Representer::ValueType>  PointIdValuePairType;
	typedef std::list<PointIdValuePairType> PointIdValueListType;
	
	typedef Domain<typename Representer::PointType> DomainType;
	
	%newobject Create;
     static StatisticalModel* Create(const Representer* representer,
									 const statismo::VectorType& m,
									 const statismo::MatrixType& orthonormalPCABasis,
									 const statismo::VectorType& pcaVariance,
									 double noiseVariance);
     ~StatisticalModel();
	 
	 %newobject Load;
     static StatisticalModel* Load(const std::string& filename, unsigned numComponents=10000);
     static StatisticalModel* Load(const H5::Group&  modelroot, unsigned numComponents=10000);
	 void Save(const std::string& filename);
	 void Save(const H5::Group& modelroot);

	const Representer* GetRepresenter() const;
	const DomainType& GetDomain() const;

	 Representer::ValueType EvaluateSampleAtPoint(Representer::DatasetConstPointerType sample, const Representer::PointType& point) const;
	DatasetPointerType DatasetToSample(DatasetConstPointerType ds) const;
	 DatasetPointerType DrawMean() const;	
	 Representer::ValueType DrawMeanAtPoint(const Representer::PointType& pt) const;
	 %rename("DrawMeanAtPointId") DrawMeanAtPoint(unsigned) const;
	 Representer::ValueType DrawMeanAtPoint(unsigned ptId) const;

	 %rename("DrawRandomSample") DrawSample() const;	 
	 DatasetPointerType DrawSample() const;	  
	 DatasetPointerType DrawSample(const statismo::VectorType& coeffs) const;
	 Representer::ValueType DrawSampleAtPoint(const statismo::VectorType& coeffs, const Representer::PointType& pt) const;	 
	 %rename("DrawSampleAtPointId") DrawSampleAtPoint(const statismo::VectorType&, unsigned) const;
	 Representer::ValueType DrawSampleAtPoint(const statismo::VectorType& coeffs, unsigned ptId) const;

   DatasetPointerType DrawPCABasisSample(unsigned componentNumber) const; 

	 
	 statismo::VectorType ComputeCoefficientsForDataset(DatasetConstPointerType ds) const;
   statismo::VectorType ComputeCoefficientsForSample(DatasetConstPointerType ds) const;
   statismo::VectorType ComputeCoefficientsForSampleVector(const statismo::VectorType& sample) const;
	 statismo::VectorType ComputeCoefficientsForPointValues(const PointValueListType&  pointValues) const;
	 statismo::VectorType ComputeCoefficientsForPointIDValues(const PointIdValueListType&  pointValues) const;
	 double ComputeLogProbabilityOfDataset(DatasetConstPointerType ds) const;
	 double ComputeProbabilityOfDataset(DatasetConstPointerType ds) const;
	 
	 statismo::MatrixType GetCovarianceAtPoint(const Representer::PointType& pt1, const Representer::PointType& pt2) const;
	 statismo::MatrixType GetCovarianceAtPoint(unsigned ptId1, unsigned ptId2) const;

	 unsigned GetNumberOfPrincipalComponents();	 	
	 statismo::VectorType DrawSampleVector(const statismo::VectorType& coefficients) const;
	 const statismo::MatrixType& GetPCABasisMatrix() const ;
	 const statismo::MatrixType GetOrthonormalPCABasisMatrix() const ;
	 const statismo::VectorType& GetPCAVarianceVector() const;
	 const statismo::VectorType& GetMeanVector() const;	 	 

	  const ModelInfo& GetModelInfo() const;
	  void SetModelInfo(const ModelInfo& modelInfo);
	double GetNoiseVariance() const;
	
	private:
	StatisticalModel();

};
} 

%template(StatisticalModel_tvr) statismo::StatisticalModel<TrivialVectorialRepresenter>;
%template(StatisticalModel_vtkPD) statismo::StatisticalModel<vtkPolyDataRepresenter>;
%template(StatisticalModel_vtkUG) statismo::StatisticalModel<vtkUnstructuredGridRepresenter>;
%template(StatisticalModel_vtkSPF3) statismo::StatisticalModel<vtkStructuredPointsRepresenter<float, 3> >;
%template(StatisticalModel_vtkSPSS1) statismo::StatisticalModel<vtkStructuredPointsRepresenter<signed short, 1> >;

//////////////////////////////////////////////////////
// PCAModelBuilder
//////////////////////////////////////////////////////


namespace statismo { 
%newobject *::BuildNewModel;
template <typename Representer>
class PCAModelBuilder {
public:

	typedef ModelBuilder<Representer> Superclass;
	typedef typename DataManager<Representer> DataManagerType;
	typedef typename DataManagerType::SampleDataStructureListType SampleDataStructureListType;	

	%newobject Create;
	static PCAModelBuilder* Create();
	
	 StatisticalModel<Representer>* BuildNewModel(const SampleDataStructureListType& sampleList, double noiseVariance, bool computeScores=true) const;	 
private:
	PCAModelBuilder();

};
}

%template(PCAModelBuilder_tvr) statismo::PCAModelBuilder<TrivialVectorialRepresenter>;
%template(PCAModelBuilder_vtkPD) statismo::PCAModelBuilder<vtkPolyDataRepresenter>;
%template(PCAModelBuilder_vtkUG) statismo::PCAModelBuilder<vtkUnstructuredGridRepresenter>;
%template(PCAModelBuilder_vtkSPF3) statismo::PCAModelBuilder<vtkStructuredPointsRepresenter<float, 3> >;
%template(PCAModelBuilder_vtkSPSS1) statismo::PCAModelBuilder<vtkStructuredPointsRepresenter<signed short, 1> >;


//////////////////////////////////////////////////////
// PartiallyFixedModelBuilder
//////////////////////////////////////////////////////



namespace statismo { 
%newobject *::BuildNewModelFromModel;
%newobject *::BuildNewModel;

template <typename Representer>
class PartiallyFixedModelBuilder {
	typedef ModelBuilder<Representer> Superclass;
public:
	typedef  DataManager<Representer> 				DataManagerType;
	typedef typename DataManagerType::SampleDataStructureListType SampleDataStructureListType;		
	typedef  StatisticalModel<Representer> 	StatisticalModelType;	
	typedef typename StatisticalModelType::PointValueType PointValueType;
	typedef  typename StatisticalModelType::PointValueListType PointValueListType;
	
	%newobject Create;
	static PartiallyFixedModelBuilder* Create();
	virtual ~PartiallyFixedModelBuilder();

	StatisticalModelType* BuildNewModelFromModel(const StatisticalModelType* model,	const PointValueListType& pointValues,  double pointValuesNoiseVariance, bool computeScores=true) const;
	StatisticalModelType* BuildNewModel(const SampleDataStructureListType& sampleList, const PointValueListType& pointValues,  double pointValuesNoiseVariance,	double noiseVariance) const;

	private:
		PartiallyFixedModelBuilder();
};
}

%template(PartiallyFixedModelBuilder_tvr) statismo::PartiallyFixedModelBuilder<TrivialVectorialRepresenter>;
%template(PartiallyFixedModelBuilder_vtkPD) statismo::PartiallyFixedModelBuilder<vtkPolyDataRepresenter>;
%template(PartiallyFixedModelBuilder_vtkUG) statismo::PartiallyFixedModelBuilder<vtkUnstructuredGridRepresenter>;
%template(PartiallyFixedModelBuilder_vtkSPF3) statismo::PartiallyFixedModelBuilder<vtkStructuredPointsRepresenter<float, 3> >;
%template(PartiallyFixedModelBuilder_vtkSPSS1) statismo::PartiallyFixedModelBuilder<vtkStructuredPointsRepresenter<signed short, 1> >;

//////////////////////////////////////////////////////
// PosteriorModelBuilder
//////////////////////////////////////////////////////


namespace statismo { 
%newobject *::BuildNewModelFromModel;
%newobject *::BuildNewModel;

template <typename Representer>
class PosteriorModelBuilder {
	typedef ModelBuilder<Representer> Superclass;
public:
	typedef  DataManager<Representer> 				DataManagerType;
	typedef typename DataManagerType::SampleDataStructureListType SampleDataStructureListType;		
	typedef  StatisticalModel<Representer> 	StatisticalModelType;	
	typedef typename StatisticalModelType::PointValueType PointValueType;
	typedef  typename StatisticalModelType::PointValueListType PointValueListType;
	
	%newobject Create;
	static PosteriorModelBuilder* Create();
	virtual ~PosteriorModelBuilder();

	StatisticalModelType* BuildNewModelFromModel(const StatisticalModelType* model,	const PointValueListType& pointValues,  double pointValuesNoiseVariance, bool computeScores=true) const;
	StatisticalModelType* BuildNewModel(const SampleDataStructureListType& sampleList, const PointValueListType& pointValues,  double pointValuesNoiseVariance,	double noiseVariance) const;

	private:
		PosteriorModelBuilder();
};
}

%template(PosteriorModelBuilder_tvr) statismo::PosteriorModelBuilder<TrivialVectorialRepresenter>;
%template(PosteriorModelBuilder_vtkPD) statismo::PosteriorModelBuilder<vtkPolyDataRepresenter>;
%template(PosteriorModelBuilder_vtkUG) statismo::PosteriorModelBuilder<vtkUnstructuredGridRepresenter>;
%template(PosteriorModelBuilder_vtkSPF3) statismo::PosteriorModelBuilder<vtkStructuredPointsRepresenter<float, 3> >;
%template(PosteriorModelBuilder_vtkSPSS1) statismo::PosteriorModelBuilder<vtkStructuredPointsRepresenter<signed short, 1> >;




//////////////////////////////////////////////////////
// ReducedVarianceModelBuilder
//////////////////////////////////////////////////////

namespace statismo { 
%newobject *::BuildNewModelFromModel;

template <typename Representer>
class ReducedVarianceModelBuilder {
	typedef ModelBuilder<Representer> Superclass;
public:		
	typedef  StatisticalModel<Representer> 	StatisticalModelType;	
	
	%newobject Create;
	static ReducedVarianceModelBuilder* Create();
	virtual ~ReducedVarianceModelBuilder();

	StatisticalModelType* BuildNewModelFromModel(const StatisticalModelType* model,	double totalVariance, bool computeScores=true) const;

	private:
		ReducedVarianceModelBuilder();
};
}

%template(ReducedVarianceModelBuilder_tvr) statismo::ReducedVarianceModelBuilder<TrivialVectorialRepresenter>;
%template(ReducedVarianceModelBuilder_vtkPD) statismo::ReducedVarianceModelBuilder<vtkPolyDataRepresenter>;
%template(ReducedVarianceModelBuilder_vtkSPF3) statismo::ReducedVarianceModelBuilder<vtkStructuredPointsRepresenter<float, 3> >;

