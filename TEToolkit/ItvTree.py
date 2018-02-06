#!/usr/bin/env python
#coding:utf-8
# Author:  Mira Han 
# Purpose: interval tree for TE indexing 
# based on the AVL tree module developed by 
# Author:  mozman (python version)
# Purpose: avl tree module (Julienne Walker's unbounded none recursive algorithm)
# source: http://eternallyconfuzzled.com/tuts/datastructures/jsw_tut_avl.aspx
# Created: 01.05.2010
# Copyright (c) 2010-2013 by Manfred Moitzi
# License: MIT License

from __future__ import absolute_import

#from .abctree import ABCTree
from array import array
from TEToolkit.Constants import TEindex_BINSIZE

__all__ = ['ItvTree']

MAXSTACK = 32


class Node(object):
    """Internal object, represents a tree node."""
    __slots__ = ['left', 'right', 'balance', 'key', 'value']

    def __init__(self, key=None, value=None):
        self.left = None
        self.right = None
        self.key = key
        self.value = value
        self.balance = 0

    def __getitem__(self, key):
        """N.__getitem__(key) <==> x[key], where key is 0 (left) or 1 (right)."""
        return self.left if key == 0 else self.right

    def __setitem__(self, key, value):
        """N.__setitem__(key, value) <==> x[key]=value, where key is 0 (left) or 1 (right)."""
        if key == 0:
            self.left = value
        else:
            self.right = value

    def free(self):
        """Remove all references."""
        self.left = None
        self.right = None
        self.key = None
        self.value = None


class ItvNode(Node, object):
    def __init__(self,start=-1,end=-1,name=-1,parent=None,left=None,right=None):
        self.__start = start
        self.__end = end 
     #   self.binstart = binstart
     #   self.binend = binend
        self.__name = name #idx in nameIDmap
        self.__namelist = {}
        #self.parent = parent
        #self.isroot = False
        self.addValue(start,end,name)
        key = self.getKey()
        value = self.getValue() 
        super(ItvNode, self).__init__(key, value)

    def getKey(self):
        bin_startID = self.__start/TEindex_BINSIZE
        if self.__start == bin_startID * TEindex_BINSIZE :
             bin_startID -= 1
        return bin_startID

    def getValue(self):
        return self.__namelist

    def setValue(self, value):
        self.__namelist = value
        return self._namelist

    def addValue(self,start,end,name) :
        if start in self.__namelist :
            self.__namelist[start].append((name,end))
        else :
            self.__namelist[start] = [(name,end)]

    def getName(self):
        return self.__name

    def overlaps(self,start,end):
        TEnamelist = []
        for s in sorted(self.__namelist.keys()) :
            if s > end :
                break
            eles = self.__namelist[s]
            for name, e in eles :
                if start <= e and end >= s :
                    #print "overlaps: ", start, " ", end, " idx=", name, " ", s, " ", e, "len: ", end-s
                    TEnamelist.append([name,min(end,e)-max(start,s)])

        return TEnamelist
        #if start < self.__end and end >= self.__start :
        #    return True
        #else :
        #    return False





def height(node):
    return node.balance if node is not None else -1


def jsw_single(root, direction):
    other_side = 1 - direction
    save = root[other_side]
    root[other_side] = save[direction]
    save[direction] = root
    rlh = height(root.left)
    rrh = height(root.right)
    slh = height(save[other_side])
    root.balance = max(rlh, rrh) + 1
    save.balance = max(slh, root.balance) + 1
    return save


def jsw_double(root, direction):
    other_side = 1 - direction
    root[other_side] = jsw_single(root[other_side], other_side)
    return jsw_single(root, direction)


class ItvTree(object):
    """
    AVLTree implements a balanced binary tree with a dict-like interface.

    see: http://en.wikipedia.org/wiki/AVL_tree

    In computer science, an AVL tree is a self-balancing binary search tree, and
    it is the first such data structure to be invented. In an AVL tree, the
    heights of the two child subtrees of any node differ by at most one;
    therefore, it is also said to be height-balanced. Lookup, insertion, and
    deletion all take O(log n) time in both the average and worst cases, where n
    is the number of nodes in the tree prior to the operation. Insertions and
    deletions may require the tree to be rebalanced by one or more tree rotations.

    The AVL tree is named after its two inventors, G.M. Adelson-Velskii and E.M.
    Landis, who published it in their 1962 paper "An algorithm for the
    organization of information."

    AVLTree() -> new empty tree.
    AVLTree(mapping) -> new tree initialized from a mapping
    AVLTree(seq) -> new tree initialized from seq [(k1, v1), (k2, v2), ... (kn, vn)]

    """
    def __init__(self, items=None):
        """T.__init__(...) initializes T; see T.__class__.__doc__ for signature"""
        self._root = None
        self._count = 0
        if items is not None:
            self.update(items)


    def _new_node(self, start=-1,end=-1,name=-1):
        """Create a new tree node."""
        self._count += 1
        return ItvNode(start,end,name)

    def insert(self, start, end, name):
        """T.insert(key, value) <==> T[key] = value, insert key, value into tree."""
        my_new_node = self._new_node(start, end, name)
        key = my_new_node.getKey()
        value = my_new_node.getValue()

        if self._root is None:
            self._root = my_new_node
        else:
            node_stack = []  # node stack
            dir_stack = array('I')  # direction stack
            done = False
            top = 0
            node = self._root
            # search for an empty link, save path
            while True:
                if key == node.key:  # update existing item
                    node.addValue(start, end, name)
                    return
                direction = 1 if key > node.key else 0
                dir_stack.append(direction)
                node_stack.append(node)
                if node[direction] is None:
                    break
                node = node[direction]

            # Insert a new node at the bottom of the tree
            node[direction] = my_new_node

            # Walk back up the search path
            top = len(node_stack) - 1
            while (top >= 0) and not done:
                direction = dir_stack[top]
                other_side = 1 - direction
                top_node = node_stack[top]
                left_height = height(top_node[direction])
                right_height = height(top_node[other_side])

                # Terminate or rebalance as necessary */
                if left_height - right_height == 0:
                    done = True
                if left_height - right_height >= 2:
                    a = top_node[direction][direction]
                    b = top_node[direction][other_side]

                    if height(a) >= height(b):
                        node_stack[top] = jsw_single(top_node, other_side)
                    else:
                        node_stack[top] = jsw_double(top_node, other_side)

                    # Fix parent
                    if top != 0:
                        node_stack[top - 1][dir_stack[top - 1]] = node_stack[top]
                    else:
                        self._root = node_stack[0]
                    done = True

                # Update balance factors
                top_node = node_stack[top]
                left_height = height(top_node[direction])
                right_height = height(top_node[other_side])

                top_node.balance = max(left_height, right_height) + 1
                top -= 1

    def remove(self, key):
        """T.remove(key) <==> del T[key], remove item <key> from tree."""
        if self._root is None:
            raise KeyError(str(key))
        else:
            node_stack = [None] * MAXSTACK  # node stack
            dir_stack = array('I', [0] * MAXSTACK)  # direction stack
            top = 0
            node = self._root

            while True:
                # Terminate if not found
                if node is None:
                    raise KeyError(str(key))
                elif node.key == key:
                    break

                # Push direction and node onto stack
                direction = 1 if key > node.key else 0
                dir_stack[top] = direction

                node_stack[top] = node
                node = node[direction]
                top += 1

            # Remove the node
            if (node.left is None) or (node.right is None):
                # Which child is not null?
                direction = 1 if node.left is None else 0

                # Fix parent
                if top != 0:
                    node_stack[top - 1][dir_stack[top - 1]] = node[direction]
                else:
                    self._root = node[direction]
                node.free()
                self._count -= 1
            else:
                # Find the inorder successor
                heir = node.right

                # Save the path
                dir_stack[top] = 1
                node_stack[top] = node
                top += 1

                while heir.left is not None:
                    dir_stack[top] = 0
                    node_stack[top] = heir
                    top += 1
                    heir = heir.left

                # Swap data
                node.key = heir.key
                node.setValue(heir.getValue())


                # Unlink successor and fix parent
                xdir = 1 if node_stack[top - 1].key == node.key else 0
                node_stack[top - 1][xdir] = heir.right
                heir.free()
                self._count -= 1

            # Walk back up the search path
            top -= 1
            while top >= 0:
                direction = dir_stack[top]
                other_side = 1 - direction
                top_node = node_stack[top]
                left_height = height(top_node[direction])
                right_height = height(top_node[other_side])
                b_max = max(left_height, right_height)

                # Update balance factors
                top_node.balance = b_max + 1

                # Terminate or rebalance as necessary
                if (left_height - right_height) == -1:
                    break
                if (left_height - right_height) <= -2:
                    a = top_node[other_side][direction]
                    b = top_node[other_side][other_side]
                    if height(a) <= height(b):
                        node_stack[top] = jsw_single(top_node, direction)
                    else:
                        node_stack[top] = jsw_double(top_node, direction)
                    # Fix parent
                    if top != 0:
                        node_stack[top - 1][dir_stack[top - 1]] = node_stack[top]
                    else:
                        self._root = node_stack[0]
                top -= 1



    #range query
    def lookup_r(self,start,end,node) :

        if node is None :
            return (None,None)

        node_start = node.getKey()
        if end < node_start :
            return self.lookup_r(start,end,node.left)

        if start > node_start :
            return self.lookup_r(start,end,node.right)

        if start == node_start and end == node_start:
            return (node,None)

        if start == node_start and end > node_start :
            return (node,self.lookup_p(end,node.right))

        if end == node_start and start < node_start :
            return (self.lookup_p(start,node.left),node)
        return (None,None)

    #point query using start point
    def lookup_p(self, start,node):
        """
        Lookup node containing data

        @param data node data object to look up
        @param parent node's parent
        @returns node and node's parent if found or None, None
        """
        if node is None :
            return None
        if start < node.getKey():
            if node.left is None:
                return None
            return self.lookup_p(start, node.left)

        elif start > node.getKey():
            if node.right is None:
                    return None
            return self.lookup_p(start, node.right)
        else:
            return node

