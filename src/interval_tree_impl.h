#ifndef __INTERVAL_TREE_IMPL__
#define __INTERVAL_TREE_IMPL__

#include "interval_tree.h"

template <class T> 
FlatIntervalTree<T>::FlatIntervalTree(void) {

}

template <class T> 
FlatIntervalTree<T>::~FlatIntervalTree(void) {

}

// disjoint_intervals[ind] will be ready to be filled after the execution of the function
template <class T> 
void FlatIntervalTree<T>::shift_right(int ind) {
	int size = disjoint_intervals.size();
	disjoint_intervals.push_back(disjoint_intervals[size-1]);
	
	for (int i = size-1; i > ind; i--) {
		disjoint_intervals[i] = disjoint_intervals[i-1];
	}
}

// returns true if and only if fresh.end > main.end
// which means overlap will continue to the next interval
template <class T> 
bool FlatIntervalTree<T>::handle_overlap(int& cur_ind, const T& fresh) {
	IntervalInfo<T>* main = &disjoint_intervals[cur_ind];
	int new_interval_ind;

	if (main->spos < fresh.start) {
		uint32_t pre_main_epos = main->epos;
		main->epos = fresh.start - 1;
		new_interval_ind = cur_ind + 1;
		
		shift_right(new_interval_ind);
		main = &disjoint_intervals[cur_ind];

		uint32_t end = (pre_main_epos < fresh.end) ? pre_main_epos : fresh.end;
		IntervalInfo<T> overlap_int(fresh.start, end, main->seg_list, fresh);
		disjoint_intervals[new_interval_ind] = overlap_int;
		
		if (pre_main_epos < fresh.end) {
			cur_ind += 2;
			return true;
		}

		else if (pre_main_epos == fresh.end)
			return false;

		else {
			new_interval_ind++;
			shift_right(new_interval_ind);
			main = &disjoint_intervals[cur_ind];
			
			IntervalInfo<T> overlap_int(fresh.end + 1, pre_main_epos, main->seg_list);
			disjoint_intervals[new_interval_ind] = overlap_int;

			return false;
		}
	}
	else {	// start <= spos
		if (main->epos < fresh.end) {
			main->push_back(fresh);
			
			cur_ind++;
			return true;
		}

		else if (main->epos == fresh.end) {
			main->push_back(fresh);
			
			return false;
		}

		else {
			uint32_t pre_main_spos = main->spos;
			main->spos = fresh.end + 1;
			new_interval_ind = cur_ind;
			shift_right(new_interval_ind);
			main = &disjoint_intervals[cur_ind];
			
			IntervalInfo<T> overlap_int(pre_main_spos, fresh.end, main->seg_list, fresh);
			disjoint_intervals[new_interval_ind] = overlap_int;

			return false;
		}
	}
}

template <class T> 
void FlatIntervalTree<T>::build(map <T, string>& sorted_list) {
	typename map <T, string>:: iterator it;

	int j = 0;
	int curr;
	for (it = sorted_list.begin(); it != sorted_list.end(); it++) {
		while (j < disjoint_intervals.size() and it->first.start > disjoint_intervals[j].epos)
			j++;

		if (j == disjoint_intervals.size()) { // not overlapping with any previous interval
			IntervalInfo<T> temp_int(it->first);
			disjoint_intervals.push_back(temp_int);
		}
		else {
			curr = j;
			bool overlap_remained = false;
			while (curr < disjoint_intervals.size()) {
				overlap_remained = handle_overlap(curr, it->first);
				if (!overlap_remained)
					break;
			}
			if (curr == disjoint_intervals.size() and overlap_remained) {
				IntervalInfo<T> temp_int(it->first);
				temp_int.spos = disjoint_intervals[curr-1].epos + 1;
				disjoint_intervals.push_back(temp_int);
			}
		}
	}
}
// assumption: target is not less than list[0]
// input interval: [, )
// return: i if target in [i-1, i)
// =>
// closest Greater than: returned index
// closest Less than or Equal: returned index - 1
template <class T>
int FlatIntervalTree<T>::search(uint32_t target) {
	int beg = 0;
	int end = disjoint_intervals.size();
	int mid;

	while (end - beg > 1) {
		mid = (beg + end) / 2;
		if (target < disjoint_intervals[mid].spos)
			end = mid;
		else
			beg = mid;
	}

	return end;
}

template <class T>
IntervalInfo<T>* FlatIntervalTree<T>::find(uint32_t pos) {
	if (pos < disjoint_intervals[0].spos)
		return NULL;

	int ind = search(pos) - 1;	//closest start less than or equal pos
	if (ind < 0 or disjoint_intervals[ind].epos < pos)
		return NULL;
	
	return &disjoint_intervals[ind];
}

template <class T>
IntervalInfo<T>* FlatIntervalTree<T>::find_ind(uint32_t pos, int& ind) {
	ind = -1;
	if (pos < disjoint_intervals[0].spos)
		return NULL;

	ind = search(pos) - 1;	//closest start less than or equal pos
	if (ind < 0 or disjoint_intervals[ind].epos < pos)
		return NULL;
	
	return &disjoint_intervals[ind];
}

template <class T>
IntervalInfo<T>* FlatIntervalTree<T>::get_node(int ind) {
	return &disjoint_intervals[ind];
}

template <class T>
void FlatIntervalTree<T>::print() {
	cerr << "size = " << disjoint_intervals.size() << endl;
	int sum = 0;
	for (int i = 0; i < disjoint_intervals.size(); i++) {
		sum += disjoint_intervals[i].epos - disjoint_intervals[i].spos + 1;
		cerr << disjoint_intervals[i].epos - disjoint_intervals[i].spos + 1 << " [" << disjoint_intervals[i].spos << ", " << disjoint_intervals[i].epos << "] = { ";
		for (int j = 0; j < disjoint_intervals[i].seg_list.size(); j++)
			cerr << disjoint_intervals[i].seg_list[j] << ", ";
		cerr << " }\n";
	}
	cout << "Total covered area: " << sum << endl;
}

#endif
