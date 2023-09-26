#include <gtest/gtest.h>

#include "CoordGraph.h"

class CoordGraphTest : public ::testing::Test {};

using namespace pulmtln;

TEST_F(CoordGraphTest, getVertices)
{
	CoordGraph g;
	g.addVertex(1);
	g.addVertex(5);

	EXPECT_EQ(2, g.getVertices().size());
}


TEST_F(CoordGraphTest, findCycles_oriented_graph)
{
	// 0->1->2
	//  <---
	CoordGraph g; 
	g.addEdge(0, 1);
	g.addEdge(1, 2);
	g.addEdge(2, 0);

	ASSERT_EQ(1, g.findCycles().size());
	EXPECT_EQ(3, g.findCycles().front().size());
}

TEST_F(CoordGraphTest, findCycles)
{
	{
		//5 ---->---- 6
		//^			  v
		//2 ---->---- 3
		//^			  v
		//1 ----<---- 0
		//^			  v	
		//4 ----<---- 7

		CoordGraph g;
		g.addEdge(4, 1);
		g.addEdge(1, 2);
		g.addEdge(2, 5);
		g.addEdge(5, 6);
		g.addEdge(6, 3);
		g.addEdge(3, 0);
		g.addEdge(0, 7);
		g.addEdge(7, 4);
		g.addEdge(2, 3);
		g.addEdge(0, 1);
		EXPECT_EQ(4, g.findCycles().size());
	}
	{
		//5 ----<---- 6
		//v			  ^
		//2 ----<---- 3
		//v			  ^
		//1 ---->---- 0
		//v			  ^	
		//4 ---->---- 7

		CoordGraph g;
		g.addEdge(1, 4);
		g.addEdge(2, 1);
		g.addEdge(5, 2);
		g.addEdge(6, 5);
		g.addEdge(3, 6);
		g.addEdge(0, 3);
		g.addEdge(7, 0);
		g.addEdge(4, 7);
		g.addEdge(3, 2);
		g.addEdge(1, 0);
		EXPECT_EQ(4, g.findCycles().size());
	}
}