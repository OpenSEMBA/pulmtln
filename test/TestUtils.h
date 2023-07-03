#pragma once

#include <gtest/gtest.h>
#include <nlohmann/json.hpp>
#include <fstream>

using json = nlohmann::json;

static std::string testDataFolder(){ return "./testData/"; }
static std::string casesFolder()   { return testDataFolder(); }
static std::string outFolder()     { return "ParaView/"; }


static std::string smbCase(const std::string& caseName)
{
	return casesFolder() + caseName + "/" + caseName + ".smb.json";
}

static std::string getCaseName()
{
	return ::testing::UnitTest::GetInstance()->current_test_info()->name();
}

static std::string getTestCaseName() 
{
	std::string suiteName{
		::testing::UnitTest::GetInstance()->current_test_info()->test_suite_name()
	};	
	return suiteName + "." + getCaseName();
}

static json readJSON(const std::string& fn)
{
	std::ifstream stream(fn);
	json j;
	stream >> j;
	return j;
}