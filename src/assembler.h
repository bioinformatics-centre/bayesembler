
/*
assembler.h - This file is part of the Bayesembler (v1.1.1)


The MIT License (MIT)

Copyright (c) 2014 Lasse Maretty and Jonas Andreas Sibbesen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/


#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include <utils.h> 
#include <fragmentLengthModel.h>
#include <queue>
#include <boost/thread.hpp> 
#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include <api/BamAux.h>
#include <api/BamWriter.h>
#include <map>
#include <tr1/unordered_set>


#include <boost/functional/hash.hpp>

using namespace std;

class Assembler {
    
    public:
        
        Assembler(OptionsContainer);
    	void runBayesemblerThreaded();

    private:

        struct ReadId {
            
            string name;
            uint hi_tag;

            bool operator == (const ReadId & ri) const { 

                return (name == ri.name and hi_tag == ri.hi_tag);
            }
        };

        struct ReadIdHasher {
            
            size_t operator()(const ReadId & ri) const {

                return ((boost::hash<string>()(ri.name) ^ (boost::hash<uint>()(ri.hi_tag) << 1)) >> 1);
            }
        };

        struct PairInfo {
            
            uint pos;
            uint insert; // e_mpos
            string cigar;
            string m_cigar;

            bool operator == (const PairInfo & pi) const { 

                return (pos == pi.pos and insert == pi.insert and cigar == pi.cigar and m_cigar == pi.m_cigar);
            }
        };

        struct PairInfoHasher {
            
            size_t operator()(const PairInfo & pi) const {

                return (((((boost::hash<uint>()(pi.pos) ^ (boost::hash<uint>()(pi.insert) << 1)) >> 1) ^ (boost::hash<string>()(pi.cigar) >> 1)) << 1) ^ (boost::hash<string>()(pi.m_cigar) >> 1));
            }
        };

        struct SingleInfo {
            
            uint pos;
            string cigar;

            bool operator == (const SingleInfo & si) const { 

                return (pos == si.pos and cigar == si.cigar);
            }
        };

        struct SingleInfoHasher {
            
            size_t operator()(const SingleInfo & si) const {

                return ((boost::hash<uint>()(si.pos) ^ (boost::hash<string>()(si.cigar) << 1)) >> 1);
            }
        };

        typedef pair<BamTools::BamAlignment*, BamTools::BamAlignment*> ReadPair;
        typedef tr1::unordered_map <ReadId, BamTools::BamAlignment*, ReadIdHasher> ReadIDs;

        typedef map <uint, tr1::unordered_map <ReadId, BamTools::BamAlignment*, ReadIdHasher> > FirstReads;
        typedef tr1::unordered_map <PairInfo, ReadPair, PairInfoHasher > ReadPairs;     

        typedef multimap <uint, BamTools::BamAlignment*> UniqueReads;


        static bool graphComparePos(GraphInfo first, GraphInfo second);
        static bool graphCompareReadCount(pair<uint,uint> first, pair<uint,uint> second);

        void generateBamIndex(string);

        double writeUniqueReads(BamTools::BamWriter *, UniqueReads *, FirstReads *);
        void addUniqueReads(UniqueReads *, ReadPairs *);
    	string generateCigarString(vector<BamTools::CigarOp> &);
        void markDuplicates(BamTools::BamAlignment &, FirstReads *, ReadPairs *);

        vector<pair<list<GraphInfo>, string> >  generateSpliceGraphs(double *);
        void graphConstructorCallback(string, list<GraphInfo> *, string, boost::mutex *, int *, int *);        
        void grapherCallback(vector<pair<list<GraphInfo>, string> > *, deque< pair<list<string>, vector<FragmentAlignment *> > *> *, boost::mutex *, boost::mutex *, int *, int *, int, bool);
		void bayesemblerCallback(deque< pair<list<string>, vector<FragmentAlignment *> > *> *, deque<AssemblyInfo *> *, boost::mutex *, boost::mutex *, int *, double, FragmentLengthModel *);
        void fragmentLengthCallback(deque< pair<list<string>, vector<FragmentAlignment *> > *> *, map<uint,uint> *, boost::mutex *, boost::mutex *, int *, uint *);	
        void outputCallback(deque<AssemblyInfo *> *, boost::mutex *, int *, int);

        double calculateSequencingProbability(string &, string &, vector<BamTools::CigarOp> &);

    	OptionsContainer options_variables;

        struct ReadInfo {
            
            uint position;
            int fragment_length;
            vector<BamTools::CigarOp> cigar_string;
            double probability;
        };


};

#endif
    
    

