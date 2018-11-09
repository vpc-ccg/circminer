#ifndef __INTERVAL_TREE_H__
#define __INTERVAL_TREE_H__

#include <stdint.h>
#include <string>
#include <vector>
#include <map>

#include "interval_info.h"
#include "common.h"

using namespace std;

template <class T> class FlatIntervalTree {
private:
	vector < IntervalInfo<T> > disjoint_intervals;
	
	void shift_right(int ind);
	bool handle_overlap(int& cur_ind, const T& fresh);

	int search(uint32_t target);

public:
	FlatIntervalTree(void);
	~FlatIntervalTree(void);

	void build(map <T, string>& sorted_list);
	IntervalInfo<T>* find(uint32_t pos);
	IntervalInfo<T>* find_ind(uint32_t pos, int& ind);

	IntervalInfo<T>* get_node(int ind);

	void print();
};

#include "interval_tree_impl.h"

#endif //__INTERVAL_TREE_H__
