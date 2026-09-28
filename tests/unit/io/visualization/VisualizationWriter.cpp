#include <filesystem>
#include <fstream>
#include <vector>
#include <gtest/gtest.h>

#include "core/Types.hpp"
#include "equation/heateq/HeatEquation.hpp"
#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/dof/DOFOrdering.hpp"
#include "io/visualization/VTUSeriesWriter.hpp"
#include "io/visualization/VisualizationWriter.hpp"
#include "mesh/ElementFamily.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/generator/BlockMesh2D.hpp"
#include "topology/TopologicalDOF.hpp"

using namespace residuum;

// --- io::visualization::VTUSeriesWriter: filename padding + .pvd manifest content ---

class VTUSeriesWriterTest : public ::testing::Test {
protected:
	const std::string dir = TEST_OUTPUT_PATH;
};

TEST_F(VTUSeriesWriterTest, AddStepPadsFilenameToFixedWidth){

	io::visualization::VTUSeriesWriter series(dir, "series_pad_test");

	const std::string f0 = series.addStep(0, 0.0);
	const std::string f10 = series.addStep(10, 0.2);
	const std::string f100 = series.addStep(100, 2.0);

	// same digit width across steps 0, 10, 100 -- the whole point of padding
	EXPECT_NE(f0.find("series_pad_test_000000.vtu"), std::string::npos);
	EXPECT_NE(f10.find("series_pad_test_000010.vtu"), std::string::npos);
	EXPECT_NE(f100.find("series_pad_test_000100.vtu"), std::string::npos);

}

TEST_F(VTUSeriesWriterTest, WritesValidPVDManifestWithRealTimeValues){

	io::visualization::VTUSeriesWriter series(dir, "series_pvd_test");

	series.addStep(0, 0.0);
	series.addStep(10, 0.2);

	std::ifstream file(std::filesystem::path(dir) / "series_pvd_test.pvd");
	ASSERT_TRUE(file.is_open());
	const std::string content(std::istreambuf_iterator<char>(file), {});

	EXPECT_NE(content.find("<VTKFile type=\"Collection\""), std::string::npos);
	EXPECT_NE(content.find("timestep=\"0\""), std::string::npos);
	EXPECT_NE(content.find("timestep=\"0.2\""), std::string::npos);
	EXPECT_NE(content.find("file=\"series_pvd_test_000000.vtu\""), std::string::npos);
	EXPECT_NE(content.find("file=\"series_pvd_test_000010.vtu\""), std::string::npos);

}

TEST_F(VTUSeriesWriterTest, PVDStaysValidMidSeries){

	io::visualization::VTUSeriesWriter series(dir, "series_midrun_test");
	series.addStep(0, 0.0);

	std::ifstream file(std::filesystem::path(dir) / "series_midrun_test.pvd");
	const std::string content(std::istreambuf_iterator<char>(file), {});

	EXPECT_NE(content.find("</Collection>"), std::string::npos);
	EXPECT_NE(content.find("</VTKFile>"), std::string::npos);

}

// --- io::visualization::VisualizationWriter ---
class VisualizationWriterTest : public ::testing::Test {
protected:

	const Real x0 = 0.0, x1 = 1.0, y0 = 0.0, y1 = 1.0;
	const Index nx = 2, ny = 2;

	static constexpr Index nsd = 2;
	static constexpr Index Px = 1, Py = 1;

	using HeatEqBundle = equation::HeatEquation<nsd, 2, mesh::ElementFamily::Quad>;

	HeatEqBundle::Basis basis{Px, Py};
	mesh::generator::BlockMesh2D gen{nx, ny, x0, x1, y0, y1, Px, Py};
	mesh::Mesh mesh;
	std::unique_ptr<topology::TopologicalDOF<HeatEqBundle::NumDOFs>> topoDOF;
	fem::boundary::EssentialBoundaryRegistry bcRegistry;
	std::vector<Real> algField;

	const std::string dir = TEST_OUTPUT_PATH;

	void SetUp() override {

		mesh = gen.generate();
		topoDOF = std::make_unique<topology::TopologicalDOF<HeatEqBundle::NumDOFs>>(mesh, fem::dof::DOFOrdering::Interleaved);
		topoDOF->buildConstraints(basis, bcRegistry);
		algField.assign(topoDOF->numFreeDOFs(), 1.0);

	}

};

TEST_F(VisualizationWriterTest, ThrowsWhenVTKFormatRequestedForTransient){

	EXPECT_THROW(io::visualization::VisualizationWriter(dir, "viz_throw_test", io::visualization::VisualizationWriter::Format::VTK, /*transient=*/true), std::runtime_error);

}

TEST_F(VisualizationWriterTest, AllowsVTKFormatForSteady){

	EXPECT_NO_THROW(io::visualization::VisualizationWriter(dir, "viz_steady_ok_test", io::visualization::VisualizationWriter::Format::VTK, /*transient=*/false));

}

TEST_F(VisualizationWriterTest, AllowsVTUFormatForTransient){

	EXPECT_NO_THROW(io::visualization::VisualizationWriter(dir, "viz_transient_ok_test", io::visualization::VisualizationWriter::Format::VTU, /*transient=*/true));

}

TEST_F(VisualizationWriterTest, VTKFormatWritesOneFixedFilenameRegardlessOfStep){

	io::visualization::VisualizationWriter viz(dir, "viz_vtk_fixedname_test", io::visualization::VisualizationWriter::Format::VTK, false);

	const std::string f0 = viz.writeField<HeatEqBundle::NumDOFs>(mesh, *topoDOF, bcRegistry, 0, 0.0, algField.data(), {"T"}, {"K"});
	const std::string f5 = viz.writeField<HeatEqBundle::NumDOFs>(mesh, *topoDOF, bcRegistry, 5, 1.0, algField.data(), {"T"}, {"K"});

	// series
	EXPECT_EQ(f0, f5);
	EXPECT_NE(f0.find("viz_vtk_fixedname_test.vtk"), std::string::npos);

	std::ifstream file(f0);
	const std::string content(std::istreambuf_iterator<char>(file), {});

	EXPECT_NE(content.find("SCALARS T[K]"), std::string::npos);

}

TEST_F(VisualizationWriterTest, VTUFormatWritesOneNumberedFilePerStepPlusPVD){

	io::visualization::VisualizationWriter viz(dir, "viz_vtu_series_test", io::visualization::VisualizationWriter::Format::VTU, true);

	const std::string f0 = viz.writeField<HeatEqBundle::NumDOFs>(mesh, *topoDOF, bcRegistry, 0, 0.0, algField.data(), {"T"}, {"K"});
	const std::string f10 = viz.writeField<HeatEqBundle::NumDOFs>(mesh, *topoDOF, bcRegistry, 10, 0.2, algField.data(), {"T"}, {"K"});

	EXPECT_NE(f0, f10);
	EXPECT_NE(f0.find("viz_vtu_series_test_000000.vtu"), std::string::npos);
	EXPECT_NE(f10.find("viz_vtu_series_test_000010.vtu"), std::string::npos);

	std::ifstream pvd(std::filesystem::path(dir) / "viz_vtu_series_test.pvd");
	ASSERT_TRUE(pvd.is_open());

	std::ifstream field(f0);
	const std::string content(std::istreambuf_iterator<char>(field), {});
	EXPECT_NE(content.find("Name=\"T\" units=\"K\""), std::string::npos);

}

TEST_F(VisualizationWriterTest, MakeCsvWriterProducesAWorkingWriter){

	io::visualization::VisualizationWriter viz(dir, "viz_csv_test", io::visualization::VisualizationWriter::Format::VTK, false);

	const auto path = std::filesystem::path(dir) / "viz_csv_test_monitor.csv";
	{
		auto csv = viz.makeCsvWriter(path.string(), {"tick", "time", "value"});
		csv.writeRow({0.0, 0.0, 42.0});
	}

	std::ifstream file(path);
	ASSERT_TRUE(file.is_open());
	const std::string content(std::istreambuf_iterator<char>(file), {});

	EXPECT_NE(content.find("tick"), std::string::npos);
	EXPECT_NE(content.find("42"), std::string::npos);

}
