#include <filesystem>
#include <iomanip>
#include <sstream>
#include <vector>
#include <gtest/gtest.h>

#include "core/Types.hpp"
#include "equation/heateq/HeatEquation.hpp"
#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/dof/DOFOrdering.hpp"
#include "io/fieldio/FieldIO.hpp"
#include "io/fieldio/TimeSeriesValueSource.hpp"
#include "mesh/ElementFamily.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/generator/BlockMesh2D.hpp"
#include "topology/TopologicalDOF.hpp"

using namespace residuum;

class TimeSeriesValueSourceTest : public ::testing::Test {
protected:

	static constexpr Index nsd = 2;
	static constexpr Index Px = 1;
	static constexpr Index Py = 1;
	static constexpr Index dofsPerNode = 1;
	static constexpr Index stepWidth = 4; // {step:04d} pattern

	using HeatEqBundle = equation::HeatEquation<nsd, 1, mesh::ElementFamily::Quad>;

	HeatEqBundle::Basis basis{Px, Py};
	mesh::generator::BlockMesh2D gen{2, 2, 0.0, 1.0, 0.0, 1.0, Px, Py};
	mesh::Mesh mesh = gen.generate();

	topology::TopologicalDOF<dofsPerNode> topoDOF{mesh, fem::dof::DOFOrdering::Interleaved};
	fem::boundary::EssentialBoundaryRegistry emptyRegistry; // no constrained DOFs

	std::filesystem::path outDir = std::filesystem::path(TEST_OUTPUT_PATH);

	void SetUp() override {
		topoDOF.buildConstraints(basis, emptyRegistry);
	}

	std::string stepFilename(Index step) const {
		std::ostringstream oss;
		oss << "series_" << std::setw(stepWidth) << std::setfill('0') << step << ".pndf";
		return (outDir / oss.str()).string();
	}

	// writes a distinct field for the given step, time = step * 0.1
	void writeStep(Index step) {

		std::vector<Real> algField(topoDOF.numFreeDOFs());
		for (Index i = 0; i < algField.size(); ++i) {
			algField[i] = static_cast<Real>(step) * 1000.0 + static_cast<Real>(i);
		}

		io::fieldio::FieldIO::writeBinary<dofsPerNode>(mesh, topoDOF, emptyRegistry, static_cast<Real>(step) * 0.1, algField.data(), stepFilename(step));

	}

};

TEST_F(TimeSeriesValueSourceTest, LoadsCorrectFramePerStep) {

	for (Index step = 0; step < 3; ++step) {
		writeStep(step);
	}

	const auto pattern = (outDir / "series_{step:04d}.pndf").string();
	io::fieldio::TimeSeriesValueSource<dofsPerNode> src(mesh, pattern);

	for (Index step = 0; step < 3; ++step) {

		for (Index nodeID = 0; nodeID < mesh.data.numNodes; ++nodeID) {

			Real out[dofsPerNode];
			src.value(nodeID, step, static_cast<Real>(step) * 0.1, out);

			const Index adof = topoDOF.toAlgebraic(topoDOF.getNodeDOF(nodeID, 0));
			const Real expected = static_cast<Real>(step) * 1000.0 + static_cast<Real>(adof);

			EXPECT_NEAR(out[0], expected, 1e-12);

		}

		EXPECT_NEAR(src.currentFrameTime(), static_cast<Real>(step) * 0.1, 1e-12);

	}

}

TEST_F(TimeSeriesValueSourceTest, MissingPlaceholderThrows) {

	const auto pattern = (outDir / "no_placeholder.pndf").string();
	io::fieldio::TimeSeriesValueSource<dofsPerNode> src(mesh, pattern);

	Real out[dofsPerNode];
	EXPECT_THROW(src.value(0, 0, 0.0, out), std::runtime_error);

}

TEST_F(TimeSeriesValueSourceTest, MismatchedTimeThrows) {

	writeStep(0); // recorded time = 0.0

	const auto pattern = (outDir / "series_{step:04d}.pndf").string();
	io::fieldio::TimeSeriesValueSource<dofsPerNode> src(mesh, pattern);

	Real out[dofsPerNode];

	// step 0 is on disk, but the caller's live solve thinks step 0 is at a different time
	EXPECT_THROW(src.value(0, 0, 5.0, out), std::runtime_error);

	// the correct time for step 0 still works
	EXPECT_NO_THROW(src.value(0, 0, 0.0, out));

}

TEST_F(TimeSeriesValueSourceTest, InterpolateModeIsStubbed) {

	writeStep(0);

	const auto pattern = (outDir / "series_{step:04d}.pndf").string();
	io::fieldio::TimeSeriesValueSource<dofsPerNode> src(mesh, pattern, io::fieldio::TimeSeriesValueSource<dofsPerNode>::Interpolation::Interpolate);

	Real out[dofsPerNode];
	EXPECT_THROW(src.value(0, 0, 0.0, out), std::runtime_error);

}
