#ifndef SPANNING_TREE_DATA_TYPES_H
#define SPANNING_TREE_DATA_TYPES_H


//an edge
struct edge 
{
  long  idx_A;      /*index of endpoint A*/
  long  idx_B;      /*index of endpoint B*/
  float dis;        /*euclidean length of edge*/
  long  peak_index; /*which peak does this edge belong to?*/
  /*Note that an index is the position in an array, while
    an id is a unique identifier for a tracer particle */
};

//a vertex
struct vertex
{
  long idx; /*index of tracer, as sorted by density*/
  long pid; /*id of group*/
};

//a linked list of vertices
struct vertex_linked_list
{
  struct vertex_linked_list *prev;
  struct vertex_linked_list *next;
  struct vertex_linked_list *head;
  struct vertex_linked_list *tail;
  long idx;
  long pid;
};


#endif /*SPANNING_TREE_DATA_TYPES_H*/