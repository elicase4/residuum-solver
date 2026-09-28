#include "mesh/Mesh.hpp"

void residuum::mesh::Mesh::print(std::ostream& os) const {

	os << "\n========================================\n";
	os << "Mesh\n";
	os << "========================================\n";

	// metadata
	os << "\nMetadata\n";
	os << "----------------------------------------\n";

	os << "parametricDim   : " << data.parametricDim << '\n';

	os << "spatialDim      : " << data.spatialDim << '\n';

	os << "numNodes        : " << data.numNodes << '\n';

	os << "numElements     : " << data.numElements << '\n';

	os << "nodesPerElement : " << data.nodesPerElement << '\n';

	os << "facesPerElement : " << data.facesPerElement << '\n';

	os << "isIGA           : " << (isIGA() ? "true" : "false") << '\n';

	// basis order
	os << "\nBasis Order\n";
	os << "----------------------------------------\n";

	os << "[ ";
	for (Index p : data.basisOrder) {
		os << p << ' ';
	}
	os << "]\n";

	// xyz
	os << "\nCoordinates (xyz)\n";
	os << "----------------------------------------\n";

	for (Index n = 0; n < data.numNodes; ++n) {
		os << "node " << std::setw(4) << n << " : [ ";

		for (Index d = 0; d < data.spatialDim; ++d) {
			os << std::setw(12) << data.xyz[n*data.spatialDim + d] << ' ';
		}

		os << "]\n";
	}

	// ien
	os << "\nConnectivity (ien)\n";
	os << "----------------------------------------\n";

	for (Index e = 0; e < data.numElements; ++e) {

		os << "elem " << std::setw(4) << e << " : [ ";

		const Index* ien = getElementNodes(e);

		for (Index a = 0; a < data.nodesPerElement; ++a) {
			os << std::setw(4) << ien[a] << ' ';
		}

		os << "]\n";
	}

	// rng
	os << "\nBoundary Tags (rng)\n";
	os << "----------------------------------------\n";

	for (Index e = 0; e < data.numElements; ++e) {
		os << "elem " << std::setw(4) << e << " : [ ";

		const Int* rng = getBoundaryTag(e);

		for (Index f = 0; f < data.facesPerElement; ++f) {
			os << std::setw(4) << rng[f] << ' ';
		}

		os << "]\n";
	}

	/*
	// extraction operator
	if (!data.C.empty()) {

		const Index extSize = extractionMatrixSize();

		os << "\nExtraction Operators (C)\n";
		os << "----------------------------------------\n";

		for (Index e = 0; e < data.numElements; ++e) {
			os << "elem " << e << '\n';

			const Real* C = getExtractionOperator(e);

			for (Index i = 0; i < extSize; ++i) {
				os << std::setw(12) << C[i] << ' ';

				// matrix formatting
				if ((i + 1) % data.nodesPerElement == 0) {
					os << '\n';
				}
			}

			os << '\n';
		}
	}
	*/

	os << "========================================\n";
}
