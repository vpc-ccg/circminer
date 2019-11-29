/*
 * Copyright (c) <2008 - 2020>, University of Washington, Simon Fraser University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, 
 * are permitted provided that the following conditions are met:
 *   
 * Redistributions of source code must retain the above copyright notice, this list
 * of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or other
 *   materials provided with the distribution.
 * - Neither the name of the <ORGANIZATION> nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Author: 
 *        Faraz Hach (fhach AT cs DOT sfu DOT ca)
 *        Iman Sarrafi (isarrafi AT cs DOT sfu DOT ca)
 */


#ifndef __HASH_TABLE__
#define __HASH_TABLE__

#include "Common.h"

typedef struct HashTable
{
	long long hv;
	int *locs; 
} HashTable;

typedef struct 
{
	GeneralIndex *list;
} IHashTable;

char			*getRefGenome();
char			*getRefGenomeName();
int				getRefGenomeOffset();
CompressedSeq	*getCmpRefGenome();
CompressedSeq	*getCmpRefGenOrigin();
int				getRefGenLength();
int				getCmpRefGenLength();
int				initLoadingCompressedGenomeMeta(char*, ContigLen**, int*);
int				initLoadingHashTableMeta(char*, ContigLen**, int*);
int				initLoadingHashTable(char*);
HashTable		*getHashTable();
GeneralIndex	*getCandidates(int hv);
unsigned char	*getAlphabetCount();
void			rewindHashTable();
int 			getChrCnt();
char 			**getChrNames();
int				*getChrLens();
int 			getMaxChrLength();
int				generateHashTableOnDisk(char*, char*);
int				generateHashTable(char*, char*);
int				checkHashTable(char*);
int				loadHashTable(double*);
int	 			loadCompressedRefGenome(double *loadTime);
void			generateCompressedGenome(char* refGen, unsigned int refGenLength, CompressedSeq* crefGen);
unsigned int	calculateHashTableSize(unsigned int *hashTable, unsigned int maxSize);
void			setHashTablePointers(GeneralIndex* fullTable, unsigned int hashTableMaxSize, unsigned int* hashTable, IHashTable* hashTablePointer);
void			finalizeLoadingCompressedGenome();
void			finalizeLoadingHashTable();

void			*calculateHashTableOnFly(int *idp);
void			*sortHashTable(int *id);

#endif
