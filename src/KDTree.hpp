// Copyright

#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include <map>
#include "Point.hpp"

using namespace std;





template <size_t N, typename ElemType>
class KDTree {
 public:
  typedef std::pair<Point<N>, ElemType> value_type;
  // constructor
  KDTree();
  // destructor
  ~KDTree();

  KDTree(const KDTree &rhs);
  KDTree &operator=(const KDTree &rhs);
  // returns dimension size of elements in tree
  size_t dimension() const;
  // returns current number of elements in tree
  size_t size() const;
  // return TRUE if tree is empty
  bool empty() const;
  // return TRUE if point is "in" tree
  bool contains(const Point<N> &pt) const;
  // find correct position of new node and inserts it in the tree
  void insert(const Point<N> &pt, const ElemType &value);
  // returns reference to value associated with point
  // if point doesnt exist in the tree, its value is added and a reference is returned 
  ElemType &operator[](const Point<N> &pt);
  // returns reference to value associated with point
  // if point doesnt exist in the tree, it throws an exception
  ElemType &at(const Point<N> &pt);
  const ElemType &at(const Point<N> &pt) const;
  // given a point and a key, find the k points nearest to the given point
  // and returns the most common value, if theres a tie the most frequent value is returned
  ElemType knn_value(const Point<N> &key, size_t k) const;


  std::vector<ElemType> knn_query(const Point<N> &key, size_t k) const;
  

 private:
  size_t dimension_;
  // number of elements in tree
  size_t size_;
  // BQP max size
  size_t BQPmax;
  // node
  struct Node {
    Point<N> key;
    ElemType val;
    size_t level;

    Node * leftNode;
    Node * rightNode;

    Node( Point<N> pt, ElemType value, size_t level, 
          Node * ln = nullptr, Node * rn = nullptr):
          key(pt), val(value), level(level), leftNode(ln), rightNode(rn){}
  };

  Node * root;

  // helper functions
  void deleteNode(Node * node);
  Node* findNode(const Point<N>& point);
  // bounded priority queue
  multimap<double, Node *> BQP;

};

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
  dimension_ = N;
  size_ = 0;
  root = nullptr;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
  deleteNode(root);
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree &rhs) {
  // TODO(me): Fill this in.
}

template <size_t N, typename ElemType>
KDTree<N, ElemType> &KDTree<N, ElemType>::operator=(const KDTree &rhs) {
  // TODO(me): Fill this in.
  return *this;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
  return dimension_;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
  return size_;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
  if (size() == 0) return true;
  return false;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N> &pt) const {
  //return findNode(pt) != nullptr;
  //Node * tmp = findNode(pt);
  //if (tmp == nullptr) return false;
  //return true;
  Node* currentNode = root;
  while (currentNode != nullptr) {
      if(currentNode->key == pt) return true;
      else if(pt[currentNode->level % N] >= currentNode->key[currentNode->level % N])
        currentNode = currentNode->rightNode;
      else
        currentNode = currentNode->leftNode;
  }
  return false;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N> &pt, const ElemType &value) {
  Node * currentNode = root;
  Node * lastNode = nullptr;
  size_t level = 0;

  while (currentNode != nullptr)
  {
    level++;
    if (pt == currentNode->key)
    {
      currentNode->val = value;
      return;
    }

    size_t keyDimension = currentNode->level % N;
    if (pt[keyDimension] < currentNode->key[keyDimension])
    {
      lastNode = currentNode;
      currentNode = currentNode->leftNode;
    }
    else if (pt[keyDimension] >= currentNode->key[keyDimension])
    {
      lastNode = currentNode;
      currentNode = currentNode->rightNode;
    }
  }
  size_++;
  Node * nNode = new Node(pt, value, level);
  if (currentNode == root)
    root = nNode;
  else if (pt[lastNode->level % N] >= lastNode->key[lastNode->level % N])
    lastNode->rightNode = nNode;
  else
    lastNode->leftNode = nNode; 
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::operator[](const Point<N> &pt) {
  Node * currentNode = root;
  Node * lastNode = nullptr;
  size_t level = 0;

  while (currentNode != nullptr)
  {
    level++;
    if (pt == currentNode->key)
      return currentNode->val;
      

    size_t keyDimension = currentNode->level % N;
    if (pt[keyDimension] < currentNode->key[keyDimension])
    {
      lastNode = currentNode;
      currentNode = currentNode->leftNode;
    }
    else if (pt[keyDimension] >= currentNode->key[keyDimension])
    {
      lastNode = currentNode;
      currentNode = currentNode->rightNode;
    }
  }
  size_++;
  Node * nNode = new Node(pt, ElemType(), level);
  if (currentNode == root){
    root = nNode;
    return nNode->val;
  } 
  else if (pt[lastNode->level % N] >= lastNode->key[lastNode->level % N])
    lastNode->rightNode = nNode;
  else
    lastNode->leftNode = nNode;
  return nNode->val; 
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) {
  Node* currentNode = root;
  while (currentNode != nullptr) {
      if(currentNode->key == pt) return currentNode->val;
      size_t keyDimension = currentNode->level % N;
      if(pt[keyDimension] >= currentNode->key[keyDimension])
        currentNode = currentNode->rightNode;
      else
        currentNode = currentNode->leftNode;
  }
  throw out_of_range("Error KDTree.at()");
}

template <size_t N, typename ElemType>
const ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) const {
  Node* currentNode = root;
  while (currentNode != nullptr) {
      if(currentNode->key == pt) return currentNode->val;
      size_t keyDimension = currentNode->level % N;
      if(pt[keyDimension] >= currentNode->key[keyDimension])
        currentNode = currentNode->rightNode;
      else
        currentNode = currentNode->leftNode;
  }
  throw out_of_range("Error KDTree.at()");
}

template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N> &key, size_t k) const {
  //multimap<double, Node*> nearestBPQ;
  //BQPmax = k;
  ElemType new_element;
  return new_element;
}

template <size_t N, typename ElemType>
std::vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N> &key,
                                                     size_t k) const {
  // TODO(me): Fill this in.
  std::vector<ElemType> values;
  return values;
}

// helper funcs

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::deleteNode(Node * node){
  if (node == nullptr) return;
  deleteNode(node->leftNode);
  deleteNode(node->rightNode);
  delete node;
}

// throws error when invoking
template <size_t N, typename ElemType>
typename KDTree<N, ElemType>::Node* 
KDTree<N, ElemType>::findNode(const Point<N>& point){
  Node * currentNode = root;
  size_t currentDim = 0;
  while (currentNode != nullptr && currentNode->key != point)
  {
    if (point[currentDim % N] < currentNode->key[currentDim % N])
      currentNode = currentNode->leftNode;
    else
      currentNode = currentNode->rightNode;
    currentDim++;
  }
  return currentNode;
}


// TODO(me): finish the implementation of the rest of the KDTree class

#endif  // SRC_KDTREE_HPP_
