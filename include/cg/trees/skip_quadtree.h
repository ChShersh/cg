#pragma once

#include <cg/trees/compressed_quadtree.h>
#include <misc/random_utils.h>

#include <random>
#include <queue>

namespace cg
{
    util::uniform_random_real<double> skip_tree_gen(0., 1.);

    template <class Scalar>
    struct skip_quadtree
    {
        Scalar lx, ly, rx, ry;
        std::vector<compressed_quadtree<Scalar>> trees;
        double threshold;

        skip_quadtree() {}

        skip_quadtree(Scalar lx, Scalar ly, Scalar rx, Scalar ry)
            : lx(lx), ly(ly), rx(rx), ry(ry), threshold(0.5)
        {
            trees.push_back(compressed_quadtree<Scalar>(lx, ly, rx, ry));
        }

        void insert(const point_2t<Scalar> & p)
        {
            Mask prev_root = "";
            std::vector<Mask> loc_positions;

            for (int i = trees.size() - 1; i >= 0; i--) {
                prev_root = trees[i].find_lowest_interesting(prev_root, p);
                loc_positions.push_back(prev_root);
            }
            trees[0].insert_from_node(loc_positions.back(), p);

            int curpos = 1;
            while (true) {
                double q = skip_tree_gen();
                if (q >= threshold) {
                    if (curpos == trees.size()) {
                        trees.push_back(compressed_quadtree<Scalar>(lx, ly, rx, ry));
                        trees[curpos].insert(p);
                        break;
                    } else {
                        trees[curpos].insert_from_node(loc_positions[loc_positions.size() - curpos - 1], p);
                        curpos++;
                    }
                } else {
                    break;
                }
            }
        }

        std::vector<std::shared_ptr<QuadNode<Scalar>>> search_all_levels(const point_2t<Scalar> & p)
        {
            Mask prev_root = "";
            std::vector<std::shared_ptr<QuadNode<Scalar>>> res(trees.size());

            for (int i = trees.size() - 1; i >= 0; i--) {
                prev_root = trees[i].find_lowest_interesting(prev_root, p);
                res[i] = trees[i].compressed_nodes[prev_root]->find(p);
            }

            return res;
        }

        std::shared_ptr<QuadNode<Scalar>> find(const point_2t<Scalar> & p)
        {
            Mask last_root = "";
            for (int i = trees.size() - 1; i >= 0; i--) {
                last_root = trees[i].find_lowest_interesting(last_root, p);
            }
            return trees[0].compressed_nodes[last_root]->find(p);
        }

        rectangle_2t<Scalar> node_rect(std::shared_ptr<QuadNode<Scalar>> node) const
        {
            return rectangle_2t<Scalar>(
                        range_t<Scalar>(node->lx, node->rx),
                        range_t<Scalar>(node->ly, node->ry)
                   );
        }

        bool is_critical(const rectangle_2t<Scalar> & eps_rect,
                         const rectangle_2t<Scalar> & quad_rect,
                         const std::shared_ptr<QuadNode<Scalar>> & node) const
        {
            for (int i = 0; i < 4; i++) {
                if (node->children[i]) {
                    auto child_rect = node_rect(node->children[i]);

                    if ((child_rect & eps_rect) == (quad_rect & eps_rect)) {
                        return false;
                    }
                }
            }
            return true;
        }

        std::shared_ptr<QuadNode<Scalar>> find_lowest_critical(const rectangle_2t<Scalar> & eps_rect,
                                                               const rectangle_2t<Scalar> & quad_rect,
                                                               const std::shared_ptr<QuadNode<Scalar>> & node,
                                                               int level)
        {
            int last_non_critical = level;
            for (int i = trees.size() - 1; i >= level + 1; i--) {
                auto qnode_it = trees[i].compressed_nodes.find(node->my_mask);

                if (qnode_it != trees[i].compressed_nodes.end() &&
                    !is_critical(eps_rect, quad_rect, qnode_it->second))
                {
                    last_non_critical = i;
                    break;
                }
            }

            auto last_node = trees[last_non_critical].compressed_nodes[node->my_mask];
            rectangle_2t<Scalar> child_rect;
            while (true) {
                bool level_back = true;
                for (int i = 0; i < 4; i++) {
                    if (last_node->children[i]) {
                        child_rect = node_rect(last_node->children[i]);

                        if ((child_rect & eps_rect) == (quad_rect & eps_rect)) {
                            if (last_non_critical == level) {
                                level_back = false;
                                last_node = last_node->children[i];
                            } else if (!(level_back = last_node->children[i]->is_leaf)) {
                                last_node = last_node->children[i];
                            }
                            break;
                        }
                    }
                }

                if (level_back) {
                    if (last_non_critical == level) break;
                    else {
                        last_non_critical--;
                        last_node = trees[last_non_critical].compressed_nodes[node->my_mask];
                    }
                }
            }

            return last_node;
        }

        void approx_rect_query(const rectangle_2t<Scalar> & rect, Scalar eps,
                               std::vector<point_2t<Scalar>> & output, int level)
        {
            if ((node_rect(trees[level].root) & rect).is_empty()) return;

            auto eps_rect = rectangle_2t<Scalar>(
                range_t<Scalar>(rect.x.inf - eps, rect.x.sup + eps),
                range_t<Scalar>(rect.y.inf - eps, rect.y.sup + eps)
            );

            std::queue<std::shared_ptr<QuadNode<Scalar>>> q;
            q.push(trees[level].root);

            while (!q.empty()) {
                auto node = q.front();
                q.pop();

                auto quad_rect = node_rect(node);

                if (node->is_leaf) {
                    if (node->point != boost::none && rect.contains(node->point.get())) {
                        output.push_back(node->point.get());
                    }
                } else if ((eps_rect & quad_rect) == quad_rect) {
                    node->add_all_subtree(output);
                } else if (!is_critical(eps_rect, quad_rect, node)) {
                    q.push(find_lowest_critical(eps_rect, quad_rect, node, level));
                } else {
                    for (int i = 0; i < 4; i++) {
                        if (node->children[i] && !(node_rect(node->children[i]) & rect).is_empty()) {
                            q.push(node->children[i]);
                        }
                    }
                }
            }
        }
    };
}
