#include <algorithm>

#include "fragment_list.h"

FragmentList::FragmentList(void) {
	head = NULL;
	tail = NULL;
	size = 0;
	max_frag_size = FRAGLIM;
}

FragmentList::FragmentList(int flen) {
	head = NULL;
	tail = NULL;
	size = 0;
	max_frag_size = FRAGLIM;
	for (int i = 0; i < flen; i++) {
		MatchedKmer* mk = new MatchedKmer();
		add_back(mk);
	}
}

FragmentList::~FragmentList(void) {
	MatchedKmer* pre_head;
	while (head != NULL) {
		pre_head = head;
		head = head->next;
		delete pre_head;
	}
}

MatchedKmer* FragmentList::get_head(void) {
	return head;
}

MatchedKmer* FragmentList::get_tail(void) {
	return tail;
}
	
int FragmentList::get_size(void) {
	return size;
}

int FragmentList::get_max_frag_size(void) {
	return max_frag_size;
}

void FragmentList::add_front(fragment_t* frags, uint32_t frag_count) {
	MatchedKmer* new_mk = new MatchedKmer(frags, frag_count, head, NULL);
	if (tail == NULL)
		tail = new_mk;
	else 
		head->prev = new_mk;
	head = new_mk;
	size++;
}

void FragmentList::add_front(MatchedKmer* new_mk) {
	new_mk->next = head;
	if (tail == NULL)
		tail = new_mk;
	else 
		head->prev = new_mk;
	head = new_mk;
	size++;
}

void FragmentList::add_back(fragment_t* frags, uint32_t frag_count) {
	MatchedKmer* new_mk = new MatchedKmer(frags, frag_count, NULL, tail);
	if (tail == NULL)
		head = new_mk;
	else
		tail->next = new_mk;
	tail = new_mk;
	size++;
}

void FragmentList::add_back(MatchedKmer* new_mk) {
	new_mk->prev = tail;
	if (tail == NULL)
		head = new_mk;
	else
		tail->next = new_mk;
	tail = new_mk;
	size++;
}

void FragmentList::sort_lists(void) {
	MatchedKmer* it;
	for (it = head; it != NULL; it = it->next) {
		std::sort(it->frags, it->frags + it->frag_count);
	}
}

void FragmentList::print() {
	vafprintf(2, stderr, "LL Size: %d\n", size);
	for (MatchedKmer* t = head; t != NULL; t = t->next)
		vafprintf(2, stderr, "Frag count: %d\n", t->frag_count);
}
