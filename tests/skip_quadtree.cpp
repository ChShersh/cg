#include <gtest/gtest.h>
#include <cg/trees/quadtree.h>
#include <cg/trees/skip_quadtree.h>

#include <misc/random_utils.h>

#include <iostream>

using namespace util;

TEST(insert_speed, naive_quadtree)
{
    using cg::quadtree;
    using cg::point_2;
    auto points = util::uniform_points(10000);
    quadtree<double> * tree = new quadtree<double>(-200, -200, 200, 200);

    for (auto pt : points) {
        tree->insert(pt);
    }

    EXPECT_TRUE(true);
}


TEST(insert_speed, compressed_quadtree)
{
    using cg::compressed_quadtree;
    using cg::point_2;
    auto points = util::uniform_points(10000);
    compressed_quadtree<double> tree(-200, -200, 200, 200);

    for (auto pt : points) {
        tree.insert(pt);
    }

    //auto mask = tree.get_max_mask();
    //std::cout << "(" << mask.first << ", " << mask.second << ")" << std::endl;

    EXPECT_TRUE(true);
}

TEST(insert_speed, skip_quadtree)
{
    using cg::skip_quadtree;
    using cg::point_2;
    auto points = util::uniform_points(10000);
    skip_quadtree<double> tree(-200, -200, 200, 200);

    for (auto pt : points) {
        tree.insert(pt);
    }

    //auto mask = tree.trees[0].get_max_mask();
    //std::cout << "(" << mask.first << ", " << mask.second << ")" << std::endl;

    EXPECT_TRUE(true);
}

TEST(find_speed, naive_quadtree)
{
    using cg::quadtree;
    using cg::point_2;
    auto points = util::uniform_points(1000);

    std::set<point_2> set_points;
    quadtree<double> * tree = new quadtree<double>(-200, -200, 200, 200);

    for (auto pt : points) {
        set_points.insert(pt);
        tree->insert(pt);
    }

    bool all_have = true;
    for (auto it = set_points.begin(); it != set_points.end(); it++) {
        auto node = tree->find(*it);
        all_have &= (node->point != boost::none) && (node->point.get() == *it);
    }

    EXPECT_TRUE(all_have);
}

TEST(find_speed, compressed_quadtree)
{
    using cg::compressed_quadtree;
    using cg::point_2;
    auto points = util::uniform_points(1000);

    std::set<point_2> set_points;
    compressed_quadtree<double> tree(-200, -200, 200, 200);

    for (auto pt : points) {
        set_points.insert(pt);
        tree.insert(pt);
    }

    bool all_have = true;
    for (auto it = set_points.begin(); it != set_points.end(); it++) {
        auto node = tree.find(*it);
        all_have &= (node->point != boost::none) && (node->point.get() == *it);
    }

    EXPECT_TRUE(all_have);
}

TEST(find_speed, skip_quadtree)
{
    using cg::skip_quadtree;
    using cg::point_2;
    auto points = util::uniform_points(1000);

    std::set<point_2> set_points;
    skip_quadtree<double> tree(-200, -200, 200, 200);

    for (auto pt : points) {
        set_points.insert(pt);
        tree.insert(pt);
    }

    bool all_have = true;
    for (auto it = set_points.begin(); it != set_points.end(); it++) {
        auto levels = tree.search_all_levels(*it);
        auto node = levels[0]->find(*it);
        all_have &= (node->point != boost::none) && (node->point.get() == *it);
    }

    EXPECT_TRUE(all_have);
}

TEST(approx_rect_speed, naive_quadtree)
{
    using cg::quadtree;
    using cg::point_2;
    using cg::rectangle_2;
    using cg::range;

    auto points = util::uniform_points(10000);
    quadtree<double> * tree = new quadtree<double>(-2000, -2000, 2000, 2000);

    rectangle_2 rect = {range{-1500, 1500}, range{-1500, 1500}};

    std::set<point_2> in_rect;
    for (auto pt : points) {
        tree->insert(pt);
        if (rect.contains(pt)) in_rect.insert(pt);
    }

    std::vector<point_2> output;
    tree->rectangle_query(rect, 10, output);
    std::set<point_2> output_set(output.begin(), output.end());

    bool all_have = true;
    for (auto pt : in_rect) {
        all_have &= (output_set.find(pt) != output_set.end());
    }

    EXPECT_TRUE(all_have);
}

TEST(approx_rect_speed, compressed_quadtree)
{
    using cg::compressed_quadtree;
    using cg::point_2;
    using cg::rectangle_2;
    using cg::range;

    auto points = util::uniform_points(10000);
    compressed_quadtree<double>  tree(-2000, -2000, 2000, 2000);

    rectangle_2 rect = {range{-1500, 1500}, range{-1500, 1500}};

    std::set<point_2> in_rect;
    for (auto pt : points) {
        tree.insert(pt);
        if (rect.contains(pt)) in_rect.insert(pt);
    }

    std::vector<point_2> output;
    tree.rectangle_query(rect, 10, output);
    std::set<point_2> output_set(output.begin(), output.end());

    bool all_have = true;
    for (auto pt : in_rect) {
        all_have &= (output_set.find(pt) != output_set.end());
    }


    EXPECT_TRUE(all_have);
}

TEST(approx_rect_speed, skip_quadtree)
{
    using cg::skip_quadtree;
    using cg::point_2;
    using cg::rectangle_2;
    using cg::range;

    auto points = util::uniform_points(10000);
    skip_quadtree<double>  tree(-2000, -2000, 2000, 2000);

    rectangle_2 rect = {range{-1500, 1500}, range{-1500, 1500}};

    std::set<point_2> in_rect;
    for (auto pt : points) {
        tree.insert(pt);
        if (rect.contains(pt)) in_rect.insert(pt);
    }

    std::vector<point_2> output;
    tree.approx_rect_query(rect, 10, output, 0);
    std::set<point_2> output_set(output.begin(), output.end());

    bool all_have = true;
    for (auto pt : in_rect) {
        all_have &= (output_set.find(pt) != output_set.end());
    }


    EXPECT_TRUE(all_have);
}
