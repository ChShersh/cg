#pragma once

#include <boost/optional.hpp>

#include <cg/primitives/point.h>
#include <cg/primitives/rectangle.h>
#include <cg/io/point.h>

namespace cg
{
    template <class Scalar>
    struct quadtree
    {
        Scalar lx, ly, rx, ry;
        bool is_leaf;
        boost::optional<point_2t<Scalar>> point;
        std::vector<quadtree *> children;

        quadtree(Scalar lx, Scalar ly, Scalar rx, Scalar ry)
          : lx(lx), ly(ly), rx(rx), ry(ry), is_leaf(true), point(boost::none)
        {
            children.resize(4);
        }

        std::vector<Scalar> coordinates_by_id(int id)
        {
            Scalar ax = lx, ay = ly, bx = rx, by = ry;
            switch (id) {
                case 0: bx = (lx + rx) / 2; by = (ly + ry) / 2; break;
                case 1: ax = (lx + rx) / 2; by = (ly + ry) / 2; break;
                case 2: ay = (ly + ry) / 2; bx = (lx + rx) / 2; break;
                case 3: ax = (lx + rx) / 2; ay = (ly + ry) / 2; break;
            }
            return std::vector<Scalar>{ax, ay, bx, by};
        }

        bool inside_id_square(int id, const point_2t<Scalar> & p)
        {
            auto c = coordinates_by_id(id);
            Scalar ax = c[0], ay = c[1], bx = c[2], by = c[3];
            return (ax <= p.x && p.x < bx && ay <= p.y && p.y < by);
        }

        void insert(const point_2t<Scalar> & p)
        {
           if (p.x < lx || rx <= p.x || p.y < ly || ry <= p.y) {
              return;
           }

           if (is_leaf) {
               if (point == boost::none) {
                   point = p;
                   return;
               } else {
                  if (p == point.get()) return;

                  is_leaf = false;
                  point_2t<Scalar> pin = point.get();
                  point = boost::none;

                  for (int i = 0; i < 4; i++) {
                      auto c = coordinates_by_id(i);
                      children[i] = new quadtree(c[0], c[1], c[2], c[3]);
                  }

                  insert(pin); // O(1) insert
               }
           }

           for (int i = 0; i < 4; i++) {
               children[i]->insert(p);
           }
        }

        bool inside_me(const point_2t<Scalar> & p)
        {
            return (lx <= p.x && p.x < rx && ly <= p.y && p.y < ry);
        }

        quadtree * find(const point_2t<Scalar> & p)
        {
            if (is_leaf) {
                return this;
            } else {
                for (int i = 0; i < 4; i++) {
                    if (children[i]->inside_me(p)) {
                        return children[i]->find(p);
                    }
                }
                return nullptr;
            }
        }

        void check_emptyness(int & has_point, bool & all_leaves,
                             boost::optional<point_2t<Scalar>> & last, quadtree * child)
        {
            if (child->is_leaf) {
                if (child->point) {
                    has_point++;
                    last = child->point;
                }
            } else {
                all_leaves = false;
            }
        }

        void remove(const point_2t<Scalar> & p)
        {
            if (p.x < lx || rx <= p.x || p.y < ly || ry <= p.y) {
                return;
            }

            if (is_leaf) {
                if (point && p == point.get()) {
                    point = boost::none;
                }
                return;
            }

            for (int i = 0; i < 4; i++) {
                if (children[i]->inside_me(p)) {
                    children[i]->remove(p);
                }
            }

            bool all_leaves = true;
            int has_point = 0;
            boost::optional<point_2t<Scalar>> last_point;
            for (int i = 0; i < 4; i++) {
                check_emptyness(has_point, all_leaves, last_point, children[i]);
            }

            if (all_leaves && has_point <= 1) {
                is_leaf = true;
                point = last_point;

                for (int i = 0; i < 4; i++) {
                    delete children[i];
                }
            }
        }

        void add_all_subtree(std::vector<point_2t<Scalar>> & output)
        {
            if (is_leaf) {
                if (point != boost::none) {
                    output.push_back(point.get());
                }
            } else {
                for (int i = 0; i < 4; i++) {
                    children[i]->add_all_subtree(output);
                }
            }
        }

        void rectangle_query(const rectangle_2t<Scalar> & rect, Scalar eps,
                             std::vector<point_2t<Scalar>> & output)
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

        ~quadtree() {
            if (!is_leaf) {
                for (int i = 0; i < 4; i++) {
                    delete children[i];
                }
            }
        }
    };
}
