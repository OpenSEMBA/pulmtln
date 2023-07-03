#include "Parser.h"

#include <filesystem>

using json = nlohmann::json;

namespace pulmtln {

const std::map<std::string, MaterialType> LABEL_TO_MATERIAL_TYPE{
	{"PEC", MaterialType::PEC}
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
		auto tag{ mat.value().at("tag").get<int>() };
		switch (type) {
		case MaterialType::PEC:
			res.pecs.push_back({ name, tag });
			break;
		default:
			throw std::runtime_error(
				"Invalid material type"
			);
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

SolverOptions Parser::readSolverOptions() const
{
	const auto& j{ json_.at("analysis") };
	
	SolverOptions res;
	setIfExists<int>(j, res.order, "order");
	setIfExists<bool>(j, res.exportParaViewSolution, "exportParaviewSolution");

	return res;
}

}