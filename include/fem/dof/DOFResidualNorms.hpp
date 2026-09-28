#ifndef RESIDUUM_FEM_DOF_DOFRESIDUALNORMS_HPP
#define RESIDUUM_FEM_DOF_DOFRESIDUALNORMS_HPP

#include <cmath>
#include <vector>

#include "core/Types.hpp"
#include "fem/dof/DOFOrdering.hpp"

namespace residuum {
	namespace fem {
		namespace dof {

			// splits a flat residual vector into a 2-norm per DOF component
			template<typename DataT>
			std::vector<DataT> computePerDOFNorms(const DataT* r, Index totalSize, Index numDOFs, Index freeDOFsPerField, DOFOrdering ordering) {

				std::vector<DataT> norms(numDOFs > 0 ? numDOFs : 1, DataT(0));

				// single field: the whole residual belongs to the one DOF
				if (numDOFs <= 1 || freeDOFsPerField == 0) {
					DataT sum = DataT(0);
					for (Index i = 0; i < totalSize; ++i) sum += r[i] * r[i];
					norms[0] = std::sqrt(sum);
					return norms;
				}

				if (ordering == DOFOrdering::Interleaved) {

					// r[i] belongs to component (i % numDOFs)
					for (Index i = 0; i < totalSize; ++i) {
						Index comp = (i % numDOFs);
						norms[comp] += r[i] * r[i];
					}

				} else {

					// r[c*freeDOFsPerField .. (c+1)*freeDOFsPerField]
					for (Index c = 0; c < numDOFs; ++c) {
						Index start = c * freeDOFsPerField;
						Index end = start + freeDOFsPerField;
						if (end > totalSize) end = totalSize;
						for (Index i = start; i < end; ++i) {
							norms[c] += r[i] * r[i];
						}
					}

				}

				for (auto& v : norms) v = std::sqrt(v);

				return norms;

			}

		} // namespace dof
	} // namespace fem
} // namespace residuum

#endif
