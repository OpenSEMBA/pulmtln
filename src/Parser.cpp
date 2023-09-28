#include "Parser.h"

#include <filesystem>

using json = nlohmann::json;

namespace pulmtln {

enum class MaterialType {
	PEC,
	Dielectric,
	OpenBoundary,
	Vacuum
};

const std::map<std::string, MaterialType> LABEL_TO_MATERIAL_TYPE{
	{"PEC", MaterialType::PEC},
	{"Dielectric", MaterialType::Dielectric},
	{"OpenBoundary", MaterialType::OpenBoundary},
	{"Vacuum", MaterialType::Vacuum}
};

json readJSON(const std::string& fn)
{
	std::ifstream stream(fn);
	if (!stream.is_open()) {
		std::runtime_error("Unable to open file: " + fn);
	}

	json j;
	stream >> j;
	return j;
}

Parser::Parser(const std::string& filename) :
	filename_{filename},
	json_{readJSON(filename)}
{}

Materials readMaterials(const json& j)
{
	Materials res;
	for (const auto& mat : j.items()) {
		auto name{ mat.key() };
		auto type{ 
			LABEL_TO_MATERIAL_TYPE.at(
				mat.value().at("type").get<std::string>()
			)
		};
		auto attribute{ mat.value().at("tag").get<int>() };
		switch (type) {
		case MaterialType::PEC:
			res.pecs.push_back({ name, attribute });
			break;
		case MaterialType::OpenBoundary:
			res.openBoundaries.push_back({ name, attribute });
			break;
		case MaterialType::Vacuum:
			res.vacuums.push_back({ name, attribute });
			break;
		case MaterialType::Dielectric:
		{
			double epsR{ mat.value().at("eps_r").get<double>() };
			res.dielectrics.push_back({ name, attribute, epsR });
			break;
		}
		default:
			throw std::runtime_error("Invalid material type");
		}
	}
	return res;
}

Model Parser::readModel() const
{
	const auto& j{ json_.at("model") };

	auto directory{
		std::filesystem::path{ filename_}.parent_path().string()
	};
	auto gmshFilename{ 
		directory + "/" + j.at("gmshFile").get<std::string>() 
	} ;

	return Model{
		mfem::Mesh::LoadFromFile(gmshFilename.c_str()),
		readMaterials(j.at("materials"))
	};
}


template <class T>
static void setIfExists(const json& j, T& entry, std::string labelToCheck)
{
	auto const it = j.find(labelToCheck);
	if (it != j.end()) {
		entry = it->get<T>();
	}
}

DriverOptions Parser::readDriverOptions() const
{
	const auto& j{ json_.at("analysis") };
	
	DriverOptions res;
	setIfExists<int>(j,  res.solverOptions.order, "order");
	setIfExists<bool>(j, res.solverOptions.printIterations, "printIterations");

	setIfExists<bool>(j, res.exportParaViewSolution, "exportParaviewSolution");
	setIfExists<bool>(j, res.exportVisItSolution, "exportVisItSolution");
	setIfExists<std::string>(j, res.exportFolder, "exportFolder");

	return res;
}

}