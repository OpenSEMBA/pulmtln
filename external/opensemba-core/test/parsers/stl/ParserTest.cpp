
#include <gtest/gtest.h>

#include "core/geometry/element/Triangle3.h"
#include "parsers/stl/Parser.h"

using namespace std;
using namespace semba;

class ParserSTLParserTest : public ::testing::Test {
protected:
    ParserSTLParserTest() {
        stlFolder_ = "./testData/";
    }

    string stlFolder_;
	string getCaseName(const string project) const {
		return stlFolder_ + project + ".stl";
	}
};

TEST_F(ParserSTLParserTest, case_nofile) 
{
    ASSERT_ANY_THROW(parsers::STL::Parser("nofile"));
}

TEST_F(ParserSTLParserTest, case_single) 
{    
    auto mesh{ 
        parsers::STL::Parser{getCaseName("single")}.readAsUnstructuredMesh()
    };

    EXPECT_EQ(3, mesh.coords().size());
    EXPECT_EQ(1, mesh.elems().getOf<geometry::Tri3>().size());
}