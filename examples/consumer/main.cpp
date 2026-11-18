// Minimal downstream consumer of the installed Residuum package. Exercises one header-only.

#include "core/Types.hpp"
#include "mesh/generator/BlockMesh2D.hpp"
#include "utils/expression/ScalarExpression.hpp"

#include <iostream>

int main() {

	residuum::mesh::generator::BlockMesh2D generator(4, 4, 0.0, 1.0, 0.0, 1.0, 1, 1);
	const residuum::mesh::Mesh mesh = generator.generate();

	std::cout << "residuum_mesh: generated " << mesh.data.numElements << " elements, " << mesh.data.numNodes << " nodes\n";

	residuum::utils::expression::ScalarExpression expr("(1 - x) * (1 - y)");

	Real coords[3] = {0.5, 0.5, 0.0};
	Real value = 0.0;
	expr(0.0, coords, &value);

	std::cout << "residuum_expression: \"(1 - x) * (1 - y)\" at (0.5, 0.5) = " << value << "\n";

	return 0;

}
