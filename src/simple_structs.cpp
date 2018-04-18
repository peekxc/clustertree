#include "simple_structs.h"

edge::edge(int from_id, double cost): from(from_id), weight(cost){
}

bool edge::operator<(const edge& other_edge) const {
  return this->weight < other_edge.weight;
}

double_edge::double_edge(int from_id, int to_id, double cost) : from(from_id), to(to_id), weight(cost) { }

bool double_edge::operator<(const double_edge& other_edge) const {
  return this->weight < other_edge.weight;
}