#ifndef SIMPLE_STRUCTS_H
#define SIMPLE_STRUCTS_H

// Structures to do priority queue
struct edge{
  int from;
  double weight;
  edge(int from_id, double cost);
  bool operator<(const edge& edge) const;
};

// Structures to do priority queue
struct double_edge{
  int from, to;
  double weight;
  double_edge(int from_id, int to_id, double cost);
  bool operator<(const double_edge& other_edge) const;
};

#endif