#include <iostream>
#include <utility>
#include <vector>

#include "benchmark/core/CSVReporter.hpp"
#include "benchmark/core/Runner.hpp"

#include "application/heateq/HeatDispatcher.hpp"
#include "application/heateq/parser/HeatConfigParser.hpp"

#include "solver/config/LinearSolverConfig.hpp"
#include "solver/config/LoggingConfig.hpp"

using namespace residuum;
using namespace residuum::application::heateq;

namespace {

	// mesh resolution -> fixture file generated from fixtures/mesh_n*.yaml via the mesh app
	struct MeshFixture {
		Index n;
		std::string file;
	};

	const std::vector<MeshFixture> meshFixtures = {
		{8, "square_n8.pmsh"},
		{16, "square_n16.pmsh"},
		{32, "square_n32.pmsh"},
		{64, "square_n64.pmsh"},
	};

	const std::vector<std::pair<std::string, solver::config::PreconditionerConfig::Type>> preconditioners = {
		{"identity", solver::config::PreconditionerConfig::Type::Identity},
		{"jacobi", solver::config::PreconditionerConfig::Type::Jacobi},
	};

} // namespace

int main() {

	config::HeatConfig baseConfig = parser::HeatConfigParser::read("config.yaml");

	// silence every I/O side effect so the timed region measures assemble+solve only
	baseConfig.logging.linear.type = solver::config::LoggerConfig::Type::None;
	baseConfig.logging.driver.type = solver::config::LoggerConfig::Type::None;
	baseConfig.output = std::nullopt;
	baseConfig.monitors.clear();

	benchmark::core::CSVReporter reporter(std::string(BENCHMARK_RESULTS_PATH) + "/HeatEq_PreconditionerComparison.csv");

	const Index warmupIters = 1;
	const Index measureIters = 3;

	for (const auto& mesh : meshFixtures) {

		for (const auto& [preconditionerName, preconditionerType] : preconditioners) {

			config::HeatConfig cfg = baseConfig;
			cfg.mesh.file = mesh.file;
			cfg.solver.linear->preconditioner.type = preconditionerType;

			bool converged = false;

			auto result = benchmark::core::run("HeatEq_PreconditionerComparison", "heateq", {{"backend", "cpu"}, {"n", std::to_string(mesh.n)}, {"preconditioner", preconditionerName}}, warmupIters, measureIters, [&]() {
				converged = HeatDispatcher::run(cfg);
			});

			if (!converged) {
				std::cerr << "PreconditionerComparison: warning -- did not converge for n=" << mesh.n << " preconditioner=" << preconditionerName << std::endl;
			}

			reporter.append(result);

			std::cout << "n=" << mesh.n << " " << preconditionerName << ": " << result.meanSeconds << "s" << std::endl;

		}

	}

	return 0;

}
