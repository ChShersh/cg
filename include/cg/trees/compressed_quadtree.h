#pragma once

#include <boost/optional.hpp>

#include <cg/primitives/point.h>
#include <cg/primitives/rectangle.h>

#include <vector>
#include <memory>
#include <unordered_map>

namespace cg
{
    typedef std::string Mask;

    template <typename Scalar>
    struct QuadNode;

    template <class Scalar>
    struct compressed_quadtree
    {
        std::shared_ptr<QuadNode<Scalar>> root;
        std::unordered_map<Mask, std::shared_ptr<QuadNode<Scalar>>> compressed_nodes;

        compressed_quadtree() {}

        compressed_quadtree(Scalar lx, Scalar ly, Scalar rx, Scalar ry)
        {
            root = std::make_shared<QuadNode<Scalar>>(lx, ly, rx, ry);
            root->my_mask = "";
            compressed_nodes[""] = root;
        }

        void insert(const point_2t<Scalar> & p) {
            root->insert(p, compressed_nodes);
        }

        void insert_from_node(const Mask & mask,
                              const point_2t<Scalar> & p)
        {
            compressed_nodes[mask]->insert(p, compressed_nodes);
        }

        Mask find_lowest_interesting(const Mask & mask,
                                     const point_2t<Scalar> & p)
        {
            return compressed_nodes[mask]->lowest_interesting(p);
        }

        std::shared_ptr<QuadNode<Scalar>> find(const point_2t<Scalar> & p) const
        {
            return root->find(p);
        }

        void remove(const point_2t<Scalar> & p)
        {
            auto node = find_lowest_interesting("", p);
        }

        void rectangle_query(const rectangle_2t<Scalar> & rect, Scalar eps,
                             std::vector<point_2t<Scalar>> & output) const
        {
            root->rectangle_query(rect, eps, output);
        }
    };

    template <typename Scalar>
    struct QuadNode : std::enable_shared_from_this<QuadNode<Scalar>>
    {
        Scalar lx, ly, rx, ry;
        bool is_leaf;
        boost::optional<point_2t<Scalar>> point;

        Mask my_mask;
        std::vector<std::shared_ptr<QuadNode>> children;

        QuadNode(Scalar lx, Scalar ly, Scalar rx, Scalar ry)
            : lx(lx), ly(ly), rx(rx), ry(ry), is_leaf(true), point(boost::none)
        {
            children.resize(4);
        }

        inline std::vector<Scalar> coordinates_by_id(Scalar plx, Scalar ply, Scalar prx,
                                                     Scalar pry, int id) const
        {
            Scalar ax = plx, ay = ply, bx = prx, by = pry;
            switch (id) {
                case 0: bx = (plx + prx) / 2; by = (ply + pry) / 2; break;
                case 1: ax = (plx + prx) / 2; by = (ply + pry) / 2; break;
                case 2: ay = (ply + pry) / 2; bx = (plx + prx) / 2; break;
                case 3: ax = (plx + prx) / 2; ay = (ply + pry) / 2; break;
            }
            return std::vector<Scalar>{ax, ay, bx, by};
        }

        inline bool inside_id_square(Scalar plx, Scalar ply, Scalar prx, Scalar pry,
                                     int id, const point_2t<Scalar> & p) const
        {
            auto c = coordinates_by_id(plx, ply, prx, pry, id);
            Scalar ax = c[0], ay = c[1], bx = c[2], by = c[3];
            return (ax <= p.x && p.x < bx && ay <= p.y && p.y < by);
        }

        inline bool inside_me(const point_2t<Scalar> & p) const
        {
            return (lx <= p.x && p.x < rx && ly <= p.y && p.y < ry);
        }

        inline int id_from_parent(Scalar plx, Scalar ply, Scalar prx, Scalar pry,
                                  const point_2t<Scalar> & p) const
        {
            for (int i = 0; i < 4; i++) {
                if (inside_id_square(plx, ply, prx, pry, i, p))
                    return i;
            }
            return -1; // should never be called
        }

        inline Mask mask_from_parent(const Mask & mask, int id)
        {
            return mask + std::to_string(id);
        }

        Mask lowest_interesting(const point_2t<Scalar> & p) const
        {
            for (int i = 0; i < 4; i++) {
                if (children[i] && children[i]->inside_me(p) && !children[i]->is_leaf) {
                    return children[i]->lowest_interesting(p);
                }
            }
            return my_mask;
        }

        std::shared_ptr<QuadNode> find(const point_2t<Scalar> & p)
        {
            for (int i = 0; i < children.size(); i++) {
                if (children[i] && children[i]->inside_me(p))
                    return children[i]->find(p);
            }
            return this->shared_from_this();
        }

        std::shared_ptr<QuadNode> insert(const point_2t<Scalar> & p,
                                         std::unordered_map<Mask, std::shared_ptr<QuadNode>> & node_map)
        {
            if (is_leaf) {
                if (point == boost::none) {
                    point = p;
                    return this->shared_from_this();
                } else {
                    if (p == point.get()) return this->shared_from_this();

                    is_leaf = false;
                    auto pin = point.get();
                    point = boost::none;

                    for (int i = 0; i < 4; i++) {
                        if (inside_id_square(lx, ly, rx, ry, i, pin)) {
                            auto c = coordinates_by_id(lx, ly, rx, ry, i);
                            children[i] = std::make_shared<QuadNode>(c[0], c[1], c[2], c[3]);
                            children[i]->point = pin;
                            children[i]->my_mask = mask_from_parent(my_mask, i);
                            node_map[children[i]->my_mask] = children[i];
                            break;
                        }
                    }
                }
            }

            for (int i = 0; i < 4; i++) {
                if (inside_id_square(lx, ly, rx, ry, i, p)) {
                    if (!children[i]) { // no such child
                        auto c = coordinates_by_id(lx, ly, rx, ry, i);
                        children[i] = std::make_shared<QuadNode>(c[0], c[1], c[2], c[3]);
                        children[i]->my_mask = mask_from_parent(my_mask, i);
                        node_map[children[i]->my_mask] = children[i];
                        children[i]->point = p;
                    } else {
                        if (children[i]->inside_me(p)) {
                            children[i] = children[i]->insert(p, node_map);
                        } else { // change nodes
                            // most hardcore part
                            auto old_child = children[i];
                            point_2t<Scalar> old_child_pt(old_child->lx, old_child->ly);
                            Scalar plx = lx, ply = ly, prx = rx, pry = ry;
                            int old_child_id, point_id;
                            Mask new_mask = my_mask;
                            do {
                                old_child_id = id_from_parent(plx, ply, prx, pry, old_child_pt),
                                point_id     = id_from_parent(plx, ply, prx, pry, p);
                                if (old_child_id != point_id) break;
                                auto c = coordinates_by_id(plx, ply, prx, pry, point_id);
                                plx = c[0]; ply = c[1]; prx = c[2]; pry = c[3];
                                new_mask = mask_from_parent(new_mask, point_id);
                            } while (true);

                            auto new_child = std::make_shared<QuadNode>(plx, ply, prx, pry);
                            new_child->children[old_child_id] = old_child;
                            new_child->is_leaf = false;
                            new_child->my_mask = new_mask;
                            node_map[new_mask] = new_child;
                            children[i] = new_child->insert(p, node_map);
                        }
                    }
                    break;
                }
            }

            int non_empty = 0, last_ind = 0;
            for (int i = 0; i < 4; i++) {
                if (children[i]) {
                    non_empty++;
                    last_ind = i;
                }
            }

            if (non_empty == 1) {
                if (my_mask != "") node_map.erase(my_mask);
                return children[last_ind];
            }

            return this->shared_from_this();
        }

        void add_all_subtree(std::vector<point_2t<Scalar>> & output) const
        {
            if (is_leaf) {
                output.push_back(point.get());
            } else {
                for (int i = 0; i < 4; i++) {
                    if (children[i]) {
                        children[i]->add_all_subtree(output);
                    }
                }
            }
        }

        void rectangle_query(const rectangle_2t<Scalar> & rect, Scalar eps,
                             std::vector<point_2t<Scalar>> & output) const
        {
            if (is_leaf) {
                if (point != boost::none && rect.contains(point.get())) {
                    output.push_back(point.get());
                }
                return;
            }

            auto eps_rect = rectangle_2t<Scalar>(
                    range_t<Scalar>(rect.x.inf - eps, rect.x.sup + eps),
                    range_t<Scalar>(rect.y.inf - eps, rect.y.sup + eps)
            );

            for (int i = 0; i < 4; i++) {
                if (!children[i]) {
                    continue;
                }
                auto c = children[i];
                auto quad_rect = rectangle_2t<Scalar>(
                    range_t<Scalar>(c->lx, c->rx),
                    range_t<Scalar>(c->ly, c->ry)
                );

                if ((eps_rect & quad_rect) == quad_rect) {
                    children[i]->add_all_subtree(output);
                } else if (!(rect & quad_rect).is_empty()) {
                    children[i]->rectangle_query(rect, eps, output);
                }
            }
        }
    };
}
