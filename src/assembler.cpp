
/*
assembler.cpp - This file is part of the Bayesembler (v1.1.1)


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


#include <assembler.h>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#include <list>
#include <vector>
#include <time.h>
#include <boost/random.hpp>
#include <set>
#include <map>
#include <gibbsSampler.h>
#include <alignmentParser.h>
#include <sequencingModel.h>
#include <allPaths.h>
#include <utils.h>
#include <valueContainer.h>
#include <gammaGenerator.h>
#include <assignmentGenerator.h>
#include <expressionGenerator.h>
#include <simplexSizeGenerator.h>
#include <boost/program_options.hpp>
#include <Eigen/Dense>
#include <gibbsOutput.h>
#include <api/BamWriter.h>
#include <limits>



Assembler::Assembler(OptionsContainer options_variables_in) {

	options_variables = options_variables_in;
}


bool Assembler::graphComparePos(GraphInfo first, GraphInfo second) {
    
    return first.left < second.left;
}


bool Assembler::graphCompareReadCount(pair<uint,uint> first, pair<uint,uint> second) {
    
    return first.second > second.second;
}


double Assembler::writeUniqueReads(BamTools::BamWriter * writer, UniqueReads * unique_reads, FirstReads * first_reads) {

    double nd_nm_pe_reads = 0;

    uint cur_nh_tag;
    uint cur_hi_tag;

    UniqueReads::iterator uit = unique_reads->begin();

    uint left_most_position;

    if (first_reads->empty()) {

        left_most_position = numeric_limits<uint>::max();
    
    } else {

        left_most_position = first_reads->begin()->first;
    }
         
    while (uit->first < left_most_position and !unique_reads->empty()) {

        cur_hi_tag = 0;
    
        assert(uit->second->GetTag("NH", cur_nh_tag));
        assert(uit->second->GetTag("HI", cur_hi_tag) or (cur_nh_tag == 1));

        if (cur_nh_tag == 1) {

            nd_nm_pe_reads += 1;
        }

        assert(writer->SaveAlignment(*(uit->second)));

        delete uit->second;
        unique_reads->erase(uit++);
    }
         
    return nd_nm_pe_reads; 
}


void Assembler::addUniqueReads(UniqueReads * unique_reads, ReadPairs * read_pairs) {

    for (ReadPairs::iterator it = read_pairs->begin(); it != read_pairs->end(); it++) {

        unique_reads->insert(unique_reads->end(), pair<uint, BamTools::BamAlignment*>(it->second.first->Position, it->second.first));
        unique_reads->insert(unique_reads->end(), pair<uint, BamTools::BamAlignment*>(it->second.second->Position, it->second.second));
    }

    read_pairs->clear();
}


string Assembler::generateCigarString(vector<BamTools::CigarOp> & cigar_data) {

    stringstream read_cigar; 

    for (vector<BamTools::CigarOp>::iterator it = cigar_data.begin(); it != cigar_data.end(); it++) {

        if ((it->Type == 'M') or (it->Type == 'D')) {

            read_cigar << it->Length << "M";

        } else if (it->Type == 'N') {

            read_cigar << it->Length << "N";
        
        } else if (!(it->Type == 'I')) {

            cerr << "ERROR: Unhandled cigar string symbol '" << it->Type << "'!" << endl;
            exit(-1);
        }
   }

   return read_cigar.str();
}


void Assembler::markDuplicates(BamTools::BamAlignment & current_alignment, FirstReads * first_reads, ReadPairs * read_pairs) {

    uint cur_hi_tag = 0;
    current_alignment.GetTag("HI", cur_hi_tag);

    ReadId ri;
    ri.name = current_alignment.Name;
    ri.hi_tag = cur_hi_tag;

    FirstReads::iterator first_reads_it = first_reads->find(current_alignment.MatePosition);

    if (first_reads_it == first_reads->end()) {

        FirstReads::iterator cur_pos_first_reads_it = first_reads->find(current_alignment.Position);

        if (cur_pos_first_reads_it == first_reads->end()) {

            ReadIDs temp_read_ids;
            assert(temp_read_ids.insert(pair<ReadId, BamTools::BamAlignment*>(ri, new BamTools::BamAlignment(current_alignment))).second);
            assert(first_reads->insert(pair<uint, ReadIDs>(current_alignment.Position, temp_read_ids)).second);

        } else {

            assert(cur_pos_first_reads_it->second.insert(pair<ReadId, BamTools::BamAlignment*>(ri, new BamTools::BamAlignment(current_alignment))).second);
        }

    } else { 
        
        ReadIDs::iterator pos_first_reads_it = first_reads_it->second.find(ri);

        if (pos_first_reads_it == first_reads_it->second.end()) {

            assert(current_alignment.Position == current_alignment.MatePosition);
            assert(first_reads_it->second.insert(pair<ReadId, BamTools::BamAlignment*>(ri, new BamTools::BamAlignment(current_alignment))).second);            
        
        } else {

            PairInfo pi;
            pi.pos = pos_first_reads_it->second->Position;
            pi.insert = abs(pos_first_reads_it->second->InsertSize);

            pi.cigar = generateCigarString(pos_first_reads_it->second->CigarData);
            pi.m_cigar = generateCigarString(current_alignment.CigarData);

            assert(current_alignment.MatePosition == pi.pos);

            ReadPairs::iterator read_pairs_it = read_pairs->find(pi);

            if (read_pairs_it == read_pairs->end()) {

                assert(read_pairs->insert(pair<PairInfo, ReadPair>(pi, ReadPair(pos_first_reads_it->second, new BamTools::BamAlignment(current_alignment)))).second);
            
            } else {

                uint cur_first_nh_tag;
                assert(pos_first_reads_it->second->GetTag("NH", cur_first_nh_tag));

                uint cur_second_nh_tag;
                assert(current_alignment.GetTag("NH", cur_second_nh_tag));

                assert(cur_first_nh_tag == cur_second_nh_tag);

                uint prev_nh_tag;
                assert(read_pairs_it->second.first->GetTag("NH", prev_nh_tag));

                if (cur_first_nh_tag == 1 and prev_nh_tag > 1) {

                    delete read_pairs_it->second.first;
                    delete read_pairs_it->second.second;
                    
                    read_pairs_it->second = ReadPair(pos_first_reads_it->second, new BamTools::BamAlignment(current_alignment));
                
                } else {

                    delete pos_first_reads_it->second;
                }
            }
            
            first_reads_it->second.erase(pos_first_reads_it);

            if (first_reads_it->second.empty()) {

                first_reads->erase(first_reads_it);
            }
        }
    }
}


void Assembler::generateBamIndex(string bam_file_name) {

    stringstream idx_system_stream;
    idx_system_stream << "nice " << options_variables.samtools_path << " index " << bam_file_name;
    
    string idx_system_string = idx_system_stream.str();
    
    if (system(NULL) == 0) {
    
        cerr << "Command processor is not avaliable" << endl;
        exit(-1);
    
    } else {
    
        assert(system(idx_system_string.c_str()) == 0);
        assert(boost::filesystem::exists(bam_file_name));
    }
}


vector<pair<list<GraphInfo>, string> > Assembler::generateSpliceGraphs(double * adjusted_fragment_count) {

	vector <string> bam_path_split;    
    boost::split(bam_path_split, options_variables.bam_file, boost::is_any_of("/"));
    
	vector <string> bam_name_split;    
    boost::split(bam_name_split, bam_path_split.back(), boost::is_any_of("."));
    
    stringstream bam_nd_pe_unstranded_stream;   
    stringstream bam_nd_pe_plus_stream;
    stringstream bam_nd_pe_minus_stream;

    for (int i = 0; i < bam_name_split.size() - 1; i++) {  

        bam_nd_pe_unstranded_stream << bam_name_split[i];
        bam_nd_pe_plus_stream << bam_name_split[i];
        bam_nd_pe_minus_stream << bam_name_split[i];
    }

    bam_nd_pe_plus_stream << "_nd_plus.bam";
    bam_nd_pe_minus_stream << "_nd_minus.bam";
    bam_nd_pe_unstranded_stream << "_nd_unstranded.bam";

    string bam_nd_pe_plus_file_name = bam_nd_pe_plus_stream.str();
    string bam_nd_pe_minus_file_name = bam_nd_pe_minus_stream.str();
    string bam_nd_pe_unstranded_file_name = bam_nd_pe_unstranded_stream.str();
         
    assert(boost::filesystem::exists(options_variables.bam_file));  
    
    // Segment bamfile into positive and negative strands if stranded
    BamTools::BamReader reader;
    assert(reader.Open(options_variables.bam_file));
    assert(reader.IsOpen());

    unsigned long total_mapped_pair_reads = 0;
    unsigned long nd_nm_mapped_paired_end_reads = 0;

    unsigned long num_reads_paired = 0;

    vector<pair<list<GraphInfo>, string> > all_graphs; 
    int graph_count = 0;

    const BamTools::SamHeader sam_header = reader.GetHeader();
    const BamTools::RefVector references = reader.GetReferenceData();

    cout << "[" << getLocalTime() << "] Removing duplicate reads" << endl;  
    
    if (options_variables.strand_specific == "unstranded") {         
        
        BamTools::BamWriter unstranded_writer;
        assert(unstranded_writer.Open(bam_nd_pe_unstranded_file_name, sam_header, references));

        BamTools::BamAlignment current_alignment;
        
        FirstReads first_reads;
        ReadPairs read_pairs;
        UniqueReads unique_reads;

        int current_position = -1;
        int cur_reference = -1;

        while (reader.GetNextAlignment(current_alignment)) {

            assert(current_alignment.IsPaired());
            assert(!current_alignment.IsFailedQC());

            if (current_alignment.IsMapped() and current_alignment.IsPrimaryAlignment() and current_alignment.IsMateMapped()) {

                total_mapped_pair_reads += 1;

                // Check that read pair is OK
                if ((current_alignment.IsReverseStrand() != current_alignment.IsMateReverseStrand()) and (current_alignment.RefID == current_alignment.MateRefID)) {

                    // Same threshold as CEM
                    if (abs(current_alignment.InsertSize) > 700000) {

                        continue;
                    }

                    if ((current_position != current_alignment.Position) or (cur_reference != current_alignment.RefID)) {

                        num_reads_paired += read_pairs.size();         

                        addUniqueReads(&unique_reads, &read_pairs);
                        nd_nm_mapped_paired_end_reads += writeUniqueReads(&unstranded_writer, &unique_reads, &first_reads);

                        if (cur_reference != current_alignment.RefID) {
        
                            cur_reference = current_alignment.RefID;

                            assert(unique_reads.empty());
                            assert(first_reads.empty());

                        } else {

                            assert(current_position < current_alignment.Position);
                        }
                    }

                    markDuplicates(current_alignment, &first_reads, &read_pairs);
                    current_position = current_alignment.Position;
                }
            }
        }

        num_reads_paired += read_pairs.size();

        addUniqueReads(&unique_reads, &read_pairs);
        nd_nm_mapped_paired_end_reads += writeUniqueReads(&unstranded_writer, &unique_reads, &first_reads);

        assert(unique_reads.empty());
        assert(first_reads.empty());

        unstranded_writer.Close();

        generateBamIndex(bam_nd_pe_unstranded_file_name);

        assert(nd_nm_mapped_paired_end_reads % 2 == 0);
        assert(total_mapped_pair_reads % 2 == 0);

        nd_nm_mapped_paired_end_reads = nd_nm_mapped_paired_end_reads / 2;
        total_mapped_pair_reads = total_mapped_pair_reads / 2;

        // Spawn threads for construction of splice-graphs, md-string and indices for each strand (unstranded_bam_file and unstranded graphs are both modified in callback function through reference)
        string unstranded_bam_file = bam_nd_pe_unstranded_file_name;
        list<GraphInfo> unstranded_graphs; 

        boost::mutex * cout_lock = new boost::mutex;
    
        cout << "[" << getLocalTime() << "] Removed duplicates from " << total_mapped_pair_reads << " mapped read pairs" << endl;   
        cout << "[" << getLocalTime() << "] Wrote " << num_reads_paired << " read pairs used for splice-graph construction" << endl;   
        cout << "\n[" << getLocalTime() << "] Spawning graph construction thread" << endl;   

        int cem_graphs = 0;
        int stats_dot_excluded = 0;

        boost::thread unstranded_thread(boost::bind(&Assembler::graphConstructorCallback, this, unstranded_bam_file, &unstranded_graphs, ".", cout_lock, &stats_dot_excluded, &cem_graphs));
        unstranded_thread.join();

        delete cout_lock;

        pair<list<GraphInfo>, string> unstranded_pair(unstranded_graphs, unstranded_bam_file);
        graph_count += unstranded_graphs.size();
        all_graphs.push_back(unstranded_pair);
        
        cout << "\n[" << getLocalTime() << "] Parsed " << cem_graphs << " splice graph(s) from cem instance file and collapsed them to " << graph_count << " assembly graph(s) (" << stats_dot_excluded << " graph(s) excluded due to inference issues resulting from unstranded data)." << endl;
              
    } else {
            
        BamTools::BamWriter plus_writer;
        assert(plus_writer.Open(bam_nd_pe_plus_file_name, sam_header, references));
        BamTools::BamWriter minus_writer;
        assert(minus_writer.Open(bam_nd_pe_minus_file_name, sam_header, references));

        pair<BamTools::BamWriter *, BamTools::BamWriter *> writer_pair;

        if (options_variables.strand_specific == "first") {

            writer_pair.first = &plus_writer;
            writer_pair.second = &minus_writer;
        
        } else {
            
            writer_pair.first = &minus_writer;
            writer_pair.second = &plus_writer;   
        }

        BamTools::BamAlignment current_alignment;

        FirstReads first_reads_1;
        FirstReads first_reads_2;

        ReadPairs read_pairs_1;
        ReadPairs read_pairs_2;

        UniqueReads unique_reads_1;
        UniqueReads unique_reads_2;
        
        int current_position_1 = -1;
        int current_position_2 = -1;
        
        int cur_reference_1 = -1;
        int cur_reference_2 = -1;

        while (reader.GetNextAlignment(current_alignment)) {

            assert(current_alignment.IsPaired());
            assert(!current_alignment.IsFailedQC());

            if (current_alignment.IsMapped() and current_alignment.IsPrimaryAlignment() and current_alignment.IsMateMapped()) {

                total_mapped_pair_reads += 1;

                // Check that read pair is OK
                if ((current_alignment.IsReverseStrand() != current_alignment.IsMateReverseStrand()) and (current_alignment.RefID == current_alignment.MateRefID)) {

                    // Same threshold as CEM
                    if (abs(current_alignment.InsertSize) > 700000) {

                        continue;
                    }

                    // Plus strand
                    if ((current_alignment.IsFirstMate() and current_alignment.IsReverseStrand() and (current_alignment.Position >= current_alignment.MatePosition)) or (!current_alignment.IsFirstMate() and !current_alignment.IsReverseStrand() and (current_alignment.Position <= current_alignment.MatePosition))) {

                        if ((current_position_1 != current_alignment.Position) or (cur_reference_1 != current_alignment.RefID)) {

                            num_reads_paired += read_pairs_1.size(); 

                            addUniqueReads(&unique_reads_1, &read_pairs_1);
                            nd_nm_mapped_paired_end_reads += writeUniqueReads(writer_pair.first, &unique_reads_1, &first_reads_1);

                            if (cur_reference_1 != current_alignment.RefID) {
            
                                cur_reference_1 = current_alignment.RefID;

                                assert(unique_reads_1.empty());
                                assert(first_reads_1.empty());

                            } else {

                                assert(current_position_1 < current_alignment.Position);
                            }
                        }

                        markDuplicates(current_alignment, &first_reads_1, &read_pairs_1); 
                        current_position_1 = current_alignment.Position;

                    // Minus strand
                    } else if ((current_alignment.IsFirstMate() and !current_alignment.IsReverseStrand() and (current_alignment.Position <= current_alignment.MatePosition)) or (!current_alignment.IsFirstMate() and current_alignment.IsReverseStrand() and (current_alignment.Position >= current_alignment.MatePosition))) {

                        if ((current_position_2 != current_alignment.Position) or (cur_reference_2 != current_alignment.RefID)) {

                            num_reads_paired += read_pairs_2.size(); 

                            addUniqueReads(&unique_reads_2, &read_pairs_2);
                            nd_nm_mapped_paired_end_reads += writeUniqueReads(writer_pair.second, &unique_reads_2, &first_reads_2);

                            if (cur_reference_2 != current_alignment.RefID) {
            
                                cur_reference_2 = current_alignment.RefID;

                                assert(unique_reads_2.empty());
                                assert(first_reads_2.empty());

                            } else {

                                assert(current_position_2 < current_alignment.Position);
                            }
                        }

                        markDuplicates(current_alignment, &first_reads_2, &read_pairs_2); 
                        current_position_2 = current_alignment.Position;           
                    }
                }
            }
        }

        num_reads_paired += read_pairs_1.size(); 

        addUniqueReads(&unique_reads_1, &read_pairs_1);
        nd_nm_mapped_paired_end_reads += writeUniqueReads(writer_pair.first, &unique_reads_1, &first_reads_1);

        assert(unique_reads_1.empty());
        assert(first_reads_1.empty()); 

        num_reads_paired += read_pairs_2.size(); 

        addUniqueReads(&unique_reads_2, &read_pairs_2);
        nd_nm_mapped_paired_end_reads += writeUniqueReads(writer_pair.second, &unique_reads_2, &first_reads_2);

        assert(unique_reads_2.empty());
        assert(first_reads_2.empty());

        plus_writer.Close();
        minus_writer.Close();

        generateBamIndex(bam_nd_pe_plus_file_name);
        generateBamIndex(bam_nd_pe_minus_file_name);

        assert(nd_nm_mapped_paired_end_reads % 2 == 0);
        assert(total_mapped_pair_reads % 2 == 0);

        nd_nm_mapped_paired_end_reads = nd_nm_mapped_paired_end_reads / 2;
        total_mapped_pair_reads = total_mapped_pair_reads / 2;

        // Spawn threads for construction of splice-graphs, md-string and indices for each strand (unstranded_bam_file and unstranded graphs are both modified in callback function through reference)
        string plus_bam_file = bam_nd_pe_plus_file_name;
        string minus_bam_file = bam_nd_pe_minus_file_name;

        list<GraphInfo> plus_graphs; 
        list<GraphInfo> minus_graphs;
        
        boost::mutex * cout_lock = new boost::mutex;

        int cem_graphs_plus = 0;
        int cem_graphs_minus = 0;
        int stats_dot_excluded_plus = 0;
        int stats_dot_excluded_minus = 0;

        cout << "[" << getLocalTime() << "] Removed duplicates from " << total_mapped_pair_reads << " mapped read pairs" << endl;   
        cout << "[" << getLocalTime() << "] Wrote " << num_reads_paired << " read pairs used for splice-graph construction" << endl;   
        cout << "\n[" << getLocalTime() << "] Spawning 2 graph construction threads" << endl;   

        boost::thread plus_thread(boost::bind(&Assembler::graphConstructorCallback, this, plus_bam_file, &plus_graphs, "+", cout_lock, &stats_dot_excluded_plus, &cem_graphs_plus));
        boost::thread minus_thread(boost::bind(&Assembler::graphConstructorCallback, this, minus_bam_file, &minus_graphs, "-", cout_lock, &stats_dot_excluded_minus, &cem_graphs_minus));
        
        plus_thread.join();
        minus_thread.join();

        assert (stats_dot_excluded_plus == 0);
        assert (stats_dot_excluded_minus == 0);

        delete cout_lock;

        pair<list<GraphInfo>, string> plus_pair(plus_graphs, plus_bam_file);
        pair<list<GraphInfo>, string> minus_pair(minus_graphs, minus_bam_file);

        graph_count += plus_graphs.size();
        graph_count += minus_graphs.size();

        all_graphs.push_back(plus_pair);
        all_graphs.push_back(minus_pair);
        
        cout << "\n[" << getLocalTime() << "] Parsed " << cem_graphs_plus + cem_graphs_minus << " splice graph(s) from cem instance files and collapsed them to " << graph_count << " assembly graph(s)" << endl;
    }

    if (graph_count == 0) {

        assert(options_variables.strand_specific == "unstranded");
        cerr << "\nWARNING: No graphs remaining after exclusion due to inference issues resulting from unstranded data" << endl;
        exit(-1);
    }   

    assert(reader.Close());
    assert (total_mapped_pair_reads >= nd_nm_mapped_paired_end_reads);

    cout << "[" << getLocalTime() << "] " << nd_nm_mapped_paired_end_reads << " unique, non-redundant read pairs being used for quantification" << endl; 

    if (options_variables.total_num_fragments == 0) {

        *adjusted_fragment_count = nd_nm_mapped_paired_end_reads;
    
    } else {

        *adjusted_fragment_count = options_variables.total_num_fragments * double(nd_nm_mapped_paired_end_reads)/total_mapped_pair_reads;
    }

    cout << "[" << getLocalTime() << "] " << *adjusted_fragment_count << " read pairs being used for FPKM normalisation\n" << endl; 

    return all_graphs;
}


void Assembler::graphConstructorCallback(string bam_nd_pe_strand_file_name, list<GraphInfo> * output_graphs, string strand, boost::mutex * cout_lock, int * stats_dot_excluded, int * cem_graphs) {    
     
    vector <string> bam_name_split;    
    boost::split(bam_name_split, bam_nd_pe_strand_file_name, boost::is_any_of("."));

    stringstream bam_nd_pe_strand_instance_prefix_stream;
    stringstream bam_nd_pe_strand_sam_stream;
    stringstream bam_nd_pe_strand_instance_stream;
    stringstream bam_nd_pe_strand_instance_log_stream;
    stringstream bam_nd_pe_strand_idx_stream;

    for (int i = 0; i < bam_name_split.size() - 1; i++) {
        
        bam_nd_pe_strand_instance_prefix_stream << bam_name_split[i];
        bam_nd_pe_strand_sam_stream << bam_name_split[i];
        bam_nd_pe_strand_instance_stream << bam_name_split[i];
        bam_nd_pe_strand_instance_log_stream << bam_name_split[i];
        bam_nd_pe_strand_idx_stream << bam_name_split[i];
    }

    bam_nd_pe_strand_sam_stream << ".sam";
    bam_nd_pe_strand_instance_stream << ".instance";
    bam_nd_pe_strand_instance_log_stream << "_instance_log.txt";
    bam_nd_pe_strand_idx_stream << ".bam.bai";

    string bam_nd_pe_strand_instance_prefix = bam_nd_pe_strand_instance_prefix_stream.str();
    string bam_nd_pe_strand_sam_name = bam_nd_pe_strand_sam_stream.str();
    string bam_nd_pe_strand_instance_file_name = bam_nd_pe_strand_instance_stream.str();
    string bam_nd_pe_strand_instance_log_file_name = bam_nd_pe_strand_instance_log_stream.str();
    string bam_nd_pe_strand_idx_file_name = bam_nd_pe_strand_idx_stream.str();
        
    cout_lock->lock();
    cout << "[" << getLocalTime() << "] Generating splice-graphs from " << bam_nd_pe_strand_file_name << " using cem" << endl;
    cout_lock->unlock();

    stringstream sam_system_stream;
    sam_system_stream << "nice " << options_variables.samtools_path << " view " << bam_nd_pe_strand_file_name << " > " << bam_nd_pe_strand_sam_name;
    string sam_system_stream_string = sam_system_stream.str();
    
    if (system(NULL) == 0) {
    
        cerr << "Command processor is not avaliable" << endl;
        exit(-1);
    
    } else {
     
        assert(system(sam_system_stream_string.c_str()) == 0);
        assert(boost::filesystem::exists(bam_nd_pe_strand_sam_name));            
    }        

    stringstream instance_system_stream;
    instance_system_stream << "nice " << options_variables.cem_processsam_path << " --no-coverage -c 1 -g 1 -u 0.05 -d " << strand << " -o " << bam_nd_pe_strand_instance_prefix << " " << bam_nd_pe_strand_sam_name << " > " << bam_nd_pe_strand_instance_log_file_name << " 2>&1";
    string instance_system_string = instance_system_stream.str();
    
    if (system(NULL) == 0) {
    
        cerr << "Command processor is not avaliable" << endl;
        exit(-1);
    
    } else {
     
        assert(system(instance_system_string.c_str()) == 0);
        assert(boost::filesystem::exists(bam_nd_pe_strand_instance_file_name));            
    }

    if (boost::filesystem::exists(bam_nd_pe_strand_sam_name)) {
    
        stringstream rm_sam_system_stream;
        rm_sam_system_stream << "rm " << bam_nd_pe_strand_sam_name;
        string rm_sam_system_string = rm_sam_system_stream.str();

        assert(system(rm_sam_system_string.c_str()) == 0);
        assert(!boost::filesystem::exists(bam_nd_pe_strand_sam_name));         
    }

    boost::regex instance_pat("^Instance\\s+(\\d+)$", boost::regex::extended);
    boost::regex boundary_pat("^Boundary\\s+([\\w|-]+)\\s+(\\d+)\\s+(\\d+)\\s+([\\+|\\.|-])$", boost::regex::extended);
    boost::regex exon_start_pat("^Segs\\s+(\\d+)$", boost::regex::extended);
    boost::regex exon_end_pat("^Refs\\s+(\\d+)$", boost::regex::extended);
    boost::regex path_start_pat("^SGTypes\\s+(\\d+)$", boost::regex::extended);
    boost::regex path_end_pat("^PETypes\\s+(\\d+)\\s+\\d+$", boost::regex::extended);
    boost::regex read_count_pat("^Reads\\s+(\\d+)$", boost::regex::extended);

    boost::match_results<const char*> matches;
    
    bool exon_line_start_pass = false;
    bool exon_line_end_pass = false;
    bool path_line_start_pass = false;
    bool path_line_end_pass = false;
    
    int graph_count = 0;
    
    int current_left_boundary = 0;
    int current_right_boundary = 0;
    int reference_exon_no = 0;
    int reference_path_no = 0;
    int real_exon_no = 0;
    int real_path_no = 0;

    uint current_read_count;
    
    string current_ref_id = "";
    string last_ref_id = "";
    string current_ref_strand = "";

    // Save graph-list separately for each chromosome to ease downstream sort
    list<list<GraphInfo> > unsorted_graphs;
    list<GraphInfo> temp_graphs;

    tr1::unordered_map <string, int> junctions;
    tr1::unordered_set <uint> excluded_exons;

    // Bookkeep vertex in- and out-degress to identify single-path graphs for fragment length estimation  
    tr1::unordered_map <int,tr1::unordered_set<int> > vertex_in_edges;
    tr1::unordered_map <int,tr1::unordered_set<int> > vertex_out_edges;
    bool single_path_graph = true;
    
    ifstream intance_file(bam_nd_pe_strand_instance_file_name.c_str());
    assert(intance_file.is_open());
    
    string instance_line;    
    stringstream dot_file_string_stream;
    
    while (intance_file.good()) {

        getline(intance_file, instance_line);
        
        if (boost::regex_match(instance_line.c_str(), matches, instance_pat)) {
        
            graph_count++;
            assert (dot_file_string_stream.str().size() == 0);
            assert (junctions.size() == 0);
            
            dot_file_string_stream << "digraph graph_" << matches[1];
            single_path_graph = true; 
        }
            
        if (boost::regex_match(instance_line.c_str(), matches, boundary_pat)) {

            current_ref_id = matches[1];
            current_ref_strand = matches[4];
                        
            if (current_ref_strand == "+") {

                dot_file_string_stream << "_plus {\n";

            } else if (current_ref_strand == "-") {

                dot_file_string_stream << "_minus {\n";

            } else {

                dot_file_string_stream << "_dot {\n";
            }

            stringstream match_stream_left(matches[2]);
            match_stream_left >> current_left_boundary;
            
            stringstream match_stream_right(matches[3]);
            match_stream_right >> current_right_boundary;
        }

        if (boost::regex_match(instance_line.c_str(), matches, read_count_pat)) {

            stringstream read_count_ss(matches[1]);
            read_count_ss >> current_read_count;
        }
            
        if (boost::regex_match(instance_line.c_str(), matches, exon_end_pat)) {    
            
            exon_line_start_pass = false;
            exon_line_end_pass = true;
            assert (reference_exon_no == real_exon_no);
        }
    
        if (exon_line_start_pass and not exon_line_end_pass) {            
            
            vector <string> exon_split;    
            boost::split(exon_split, instance_line, boost::is_any_of("\t"));

            if (options_variables.exon_base_coverage <= (1-atof(exon_split[7].c_str()))) {
            
                int left_exon_boundary = atoi(exon_split[0].c_str());
                
                if (left_exon_boundary < current_left_boundary) {
                    
                    current_left_boundary = left_exon_boundary;
                }
                
                int right_exon_boundary = atoi(exon_split[1].c_str());
                          
                if (right_exon_boundary > current_right_boundary) {
                    
                    current_right_boundary = right_exon_boundary;
                }
                
                dot_file_string_stream << real_exon_no << " [label=\"" << current_ref_id << ";" << current_ref_strand << ";" << exon_split[0] << ";" << exon_split[1] << "\"]\n";
            
            } else {

                assert(excluded_exons.insert(real_exon_no).second);
            }

            real_exon_no += 1;    
        }

        if (boost::regex_match(instance_line.c_str(), matches, exon_start_pat)) {
            
            exon_line_start_pass = true;
            exon_line_end_pass = false;
    
            real_exon_no = 0;
            
            stringstream match_stream_exon(matches[1]);
            match_stream_exon >> reference_exon_no;                        
        }

        if (boost::regex_match(instance_line.c_str(), matches, path_end_pat)) {
            
            path_line_start_pass = false;
            path_line_end_pass = true;
            assert (reference_path_no == real_path_no);
            
            if (!options_variables.no_pre_mrna) {
                
                dot_file_string_stream << "pre_mrna [label=\"" << current_ref_id << ";" << current_ref_strand << ";" << current_left_boundary << ";" << current_right_boundary << "\"]\n";
            }
            
            for (tr1::unordered_map <string, int>::iterator it = junctions.begin(); it != junctions.end(); it++) {
                
                dot_file_string_stream << it->first << " [coverage=\"" << it->second << "\"]\n";
            }

            
            dot_file_string_stream << "}\n";
            
            if (current_ref_id != last_ref_id) {

                if (temp_graphs.size() > 0) {

                    unsorted_graphs.push_back(temp_graphs);
                    temp_graphs.clear();
                }

                last_ref_id = current_ref_id;
            } 

            if (current_ref_strand == "+" or current_ref_strand == "-") {            

                string dot_string = dot_file_string_stream.str();


                list<string> dot_strings; 
                assert(dot_string != "");
                dot_strings.push_back(dot_string);
                
                GraphInfo current_graph;
                current_graph.reference = current_ref_id;
                current_graph.strand = current_ref_strand;
                current_graph.left = current_left_boundary;
                current_graph.right = current_right_boundary;
                current_graph.single_path_graph = single_path_graph;
                current_graph.dot_strings = dot_strings; 
                current_graph.read_count = current_read_count;
                temp_graphs.push_back(current_graph);

            } else {

                assert(current_ref_strand == ".");
                assert(strand == ".");
                (*stats_dot_excluded)++;
            }

            dot_file_string_stream.clear();
            dot_file_string_stream.str(std::string());
            
            junctions.clear();
            excluded_exons.clear();
            vertex_in_edges.clear();
            vertex_out_edges.clear();
        }
                           
        if (path_line_start_pass and not path_line_end_pass) {
            
            vector <string> exon_split;    
            boost::split(exon_split, instance_line, boost::is_any_of(" "));
            
            int first_idx = -1;
            assert (exon_split.size() == real_exon_no + 1);
            
            for (int i = 0; i < real_exon_no; i++) {

                if (excluded_exons.count(i) > 0) {

                    continue;
                }
                
                if (exon_split[i][0] == '1') {
                    
                    assert (exon_split[i].size() == 1);
                    first_idx = i;
                    break;
                }
            }
            
            int last_idx = first_idx;
                       
            for (int i = first_idx + 1; i < real_exon_no; i++) {

                if (excluded_exons.count(i) > 0) {

                    continue;
                }
                                
                if (exon_split[i][0] == '1') {
                
                    assert (exon_split[i].size() == 1);
                    stringstream cur_junction_stream;        
                    cur_junction_stream << last_idx << "->" << i;
                    string cur_junction = cur_junction_stream.str();

                    // Bookkeep vertex in- and out-degress to identify single-path graphs for fragment length estimation  
                    if (single_path_graph) {

                        pair<tr1::unordered_set<int>::iterator, bool> vertex_in_insert = vertex_in_edges[i].insert(last_idx);
                        pair<tr1::unordered_set<int>::iterator, bool> vertex_out_insert = vertex_out_edges[last_idx].insert(i);

                        // Only check size if insertion was successfull
                        if (vertex_in_insert.second) {

                            if (vertex_in_edges[i].size() > 1) {
                                
                                single_path_graph = false; 
                            }
                        } 

                        // Only check size if insertion was successfull
                        if (vertex_out_insert.second) {

                            if (vertex_out_edges[last_idx].size() > 1) {
                                
                                single_path_graph = false; 
                            }
                        } 
                    }

                    last_idx = i;
                        
                    if (junctions.count(cur_junction) > 0) {
                        
                        junctions[cur_junction] += atoi(exon_split[real_exon_no].c_str());
                        
                    } else {
                        
                        junctions[cur_junction] = atoi(exon_split[real_exon_no].c_str());
                    }
                }
            }
                                
            real_path_no += 1;  
        }     
            
        if (boost::regex_match(instance_line.c_str(), matches, path_start_pat)) {
            
            path_line_start_pass = true;
            path_line_end_pass = false;
            
            real_path_no = 0;
            
            stringstream match_stream_path(matches[1]);
            match_stream_path >> reference_path_no; 
        }
    }

    if (current_ref_strand == "+" or current_ref_strand == "-") {

        assert(temp_graphs.size() > 0);

    }

    if (temp_graphs.size() > 0) {

        unsorted_graphs.push_back(temp_graphs);

    }
    
    intance_file.close();
    *cem_graphs = graph_count;

    if (strand != ".") {
        
        assert ((*stats_dot_excluded) == 0); 
    }

    cout_lock->lock();
    cout << "[" << getLocalTime() << "] Parsed " << graph_count << " graph(s) from cem instance file" << endl;
    cout_lock->unlock();

    if (!options_variables.no_pre_mrna) {

        GraphInfo last_graph_plus; 
        GraphInfo last_graph_minus;

        bool first_graph_plus = true;
        bool first_graph_minus = true;

        int collapsed_plus = 0;
        int collapsed_minus = 0;

        for (list<list<GraphInfo> >::iterator reference_iter = unsorted_graphs.begin(); reference_iter != unsorted_graphs.end(); reference_iter++) {

            reference_iter->sort(Assembler::graphComparePos);

            for (list<GraphInfo>::iterator graph_iter = reference_iter->begin(); graph_iter != reference_iter->end(); graph_iter++) {

                if (graph_iter->strand == "+") {

                    if (!first_graph_plus and graph_iter->reference == last_graph_plus.reference) {
                            
                        assert(last_graph_plus.left <= graph_iter->left);
                    }

                    // Check if graph should be merged
                    if (!first_graph_plus and graph_iter->reference == last_graph_plus.reference and graph_iter->left <= last_graph_plus.right) {

                        last_graph_plus.right = max(graph_iter->right, last_graph_plus.right);

                        // Merge dot-strings
                        assert(graph_iter->dot_strings.size() == 1);
                        last_graph_plus.single_path_graph = false;
                        last_graph_plus.dot_strings.push_back(graph_iter->dot_strings.front());
                        last_graph_plus.read_count += graph_iter->read_count;
                        collapsed_plus++;
                    
                    } else {

                        if (!first_graph_plus) {
                            
                            assert(last_graph_plus.dot_strings.size() > 0);
                            output_graphs->push_back(last_graph_plus);

                        } else {

                            first_graph_plus = false;
                        }

                        last_graph_plus = *graph_iter;
                    }

                } else {

                    assert(graph_iter->strand == "-");

                    if (!first_graph_minus and graph_iter->reference == last_graph_minus.reference) {
                            
                        assert(last_graph_minus.left <= graph_iter->left);
                    }

                    // Check if graph should be merged
                    if (!first_graph_minus and graph_iter->reference == last_graph_minus.reference and graph_iter->left <= last_graph_minus.right) {

                        last_graph_minus.right = max(graph_iter->right, last_graph_minus.right);
                        last_graph_minus.single_path_graph = false;

                        // Merge dot-strings
                        assert(graph_iter->dot_strings.size() == 1);
                        last_graph_minus.dot_strings.push_back(graph_iter->dot_strings.front());
                        last_graph_minus.read_count += graph_iter->read_count;
                        collapsed_minus++;
                    
                    } else {

                        if (!first_graph_minus) {

                            assert(last_graph_minus.dot_strings.size() > 0);
                            output_graphs->push_back(last_graph_minus);

                        } else {

                            first_graph_minus = false;
                        }

                        last_graph_minus = *graph_iter;
                    }
                }
            }
        }

        if (!first_graph_plus) {

            output_graphs->push_back(last_graph_plus);
        }

        if (!first_graph_minus) {

            output_graphs->push_back(last_graph_minus);
        }

        assert(graph_count == output_graphs->size() + (*stats_dot_excluded) + collapsed_plus + collapsed_minus);
    
    } else {

        for (list<list<GraphInfo> >::iterator reference_iter = unsorted_graphs.begin(); reference_iter != unsorted_graphs.end(); reference_iter++) {

            for (list<GraphInfo>::iterator graph_iter = reference_iter->begin(); graph_iter != reference_iter->end(); graph_iter++) {

                output_graphs->push_back(*graph_iter);
            }
        }

        assert(graph_count == output_graphs->size() + (*stats_dot_excluded));
    }
    
    if (!options_variables.keep_temp_files) {
    
        if (system(NULL) == 0) {
    
            cerr << "Command processor is not avaliable" << endl;
            exit(-1);
    
        } else {
        
            stringstream rm_instance_system_stream;
            rm_instance_system_stream << "rm " << bam_nd_pe_strand_instance_file_name;
            string rm_instance_system_string = rm_instance_system_stream.str();

            assert(system(rm_instance_system_string.c_str()) == 0);
            assert(!boost::filesystem::exists(bam_nd_pe_strand_instance_file_name));  
    
            stringstream rm_instance_log_system_stream;
            rm_instance_log_system_stream << "rm " << bam_nd_pe_strand_instance_log_file_name;
            string rm_instance_log_system_string = rm_instance_log_system_stream.str();

            assert(system(rm_instance_log_system_string.c_str()) == 0);
            assert(!boost::filesystem::exists(bam_nd_pe_strand_instance_log_file_name));  
        }
    }
}


void Assembler::grapherCallback(vector<pair<list<GraphInfo>, string> > * graph_list, deque< pair<list<string>, vector<FragmentAlignment *> > *> * graph_queue, boost::mutex * graph_lock, boost::mutex * output_lock, int * remaining_bayesembler_graphs, int * remaining_output_graphs, int max_buffer_size, bool keep_temp_files) {
    
    int graph_counter = 0;
    tr1::unordered_set <string> unique_files;
    
    vector<FragmentAlignment *> read_list;
    read_list.reserve(10000);

    tr1::unordered_map<string, ReadInfo> single_mates;

    double probability = 1;

    for (vector<pair<list<GraphInfo>, string> >::iterator it = graph_list->begin(); it != graph_list->end(); it++) {

        BamTools::BamReader region_reader;
    	assert(region_reader.Open(it->second));
        assert(region_reader.IsOpen());
        assert(region_reader.LocateIndex());      
        assert(region_reader.HasIndex()); 

        unique_files.insert(it->second);
        
        list<GraphInfo>::iterator graph_it = it->first.begin();
        
        while (graph_it != it->first.end()) {

            graph_counter += 1;
            read_list.clear();
            single_mates.clear();
                    
    	    int ref_id = region_reader.GetReferenceID(graph_it->reference);
            assert(ref_id > -1);
    	    assert(region_reader.SetRegion(ref_id, graph_it->left, ref_id, graph_it->right));
            
        	BamTools::BamAlignment current_alignment;   

            while (region_reader.GetNextAlignment(current_alignment)) {

                if (!current_alignment.IsMateMapped()) {
                    
                    continue;
                }

                uint nh_tag;
                assert(current_alignment.GetTag("NH", nh_tag));

                if (nh_tag > 1) {

                    continue;
                }

    	        if ((current_alignment.MatePosition + 1 > graph_it->right) or (current_alignment.MatePosition + 1 < graph_it->left)) {
                
    	            continue; 
    	        }

                string md_tag;
                assert(current_alignment.GetTag("MD", md_tag));

                if (single_mates.count(current_alignment.Name) > 0) {

                    ReadInfo * mate = &(single_mates[current_alignment.Name]);

                    assert (mate->position <= current_alignment.Position);
                    assert (mate->position == current_alignment.MatePosition);

                    if (mate->fragment_length != -current_alignment.InsertSize) {

                        assert(mate->position == current_alignment.Position);
                    }

                    FragmentAlignment * cur_read = new FragmentAlignment;
                    cur_read->read_data.first = ReadData(mate->position, mate->cigar_string);
                    cur_read->read_data.second = ReadData(current_alignment.Position, current_alignment.CigarData);
                    cur_read->probability = mate->probability * calculateSequencingProbability(md_tag, current_alignment.Qualities, current_alignment.CigarData);
                                                
                    read_list.push_back(cur_read);
                    assert(single_mates.erase(current_alignment.Name));

                } else {

                    ReadInfo ri;
                    ri.position = current_alignment.Position;
                    ri.fragment_length = current_alignment.InsertSize;
                    ri.cigar_string = current_alignment.CigarData;
                    ri.probability = calculateSequencingProbability(md_tag, current_alignment.Qualities, current_alignment.CigarData);

                    assert (single_mates.insert(pair<string, ReadInfo>(current_alignment.Name, ri)).second);
                }         
            }    

            if (read_list.size() > 0) {

                pair<list<string>, vector<FragmentAlignment *> > * graph_pair = new pair<list<string>, vector<FragmentAlignment *> >(graph_it->dot_strings, read_list);

                // Place graph in queue
                graph_lock->lock();
                graph_queue->push_back(graph_pair);
                graph_lock->unlock();

            } else {

                graph_lock->lock();
                (*remaining_bayesembler_graphs)--;
                graph_lock->unlock();
                
                output_lock->lock();
                (*remaining_output_graphs)--;
                output_lock->unlock();
            }
            
            graph_it = it->first.erase(graph_it);
            
        	while (graph_queue->size() > max_buffer_size) {

        		boost::this_thread::sleep(boost::posix_time::seconds(5));
        	}
        }
        
        assert(region_reader.Close());
    }

    if (!keep_temp_files) {   

        for (tr1::unordered_set<string>::iterator unique_file_it = unique_files.begin(); unique_file_it != unique_files.end(); unique_file_it++) {

            if (system(NULL) == 0) {
    
                cerr << "Command processor is not avaliable" << endl;
                exit(-1);
    
            } else {
            
                stringstream rm_bam_nd_pe_system_stream;
                rm_bam_nd_pe_system_stream << "rm " << *unique_file_it;
                string rm_bam_nd_pe_system_string = rm_bam_nd_pe_system_stream.str();

                assert(system(rm_bam_nd_pe_system_string.c_str()) == 0);
                assert(!boost::filesystem::exists(*unique_file_it));  
                
                stringstream rm_bam_nd_pe_idx_name;
                rm_bam_nd_pe_idx_name << *unique_file_it << ".bai";
                
                stringstream rm_bam_nd_pe_idx_system_stream;
                rm_bam_nd_pe_idx_system_stream << "rm " << rm_bam_nd_pe_idx_name.str();
                string rm_bam_nd_pe_idx_system_string = rm_bam_nd_pe_idx_system_stream.str();

                assert(system(rm_bam_nd_pe_idx_system_string.c_str()) == 0);
                assert(!boost::filesystem::exists(rm_bam_nd_pe_idx_name.str()));      
            }
        }
    }
}


double Assembler::calculateSequencingProbability(string & md_tag, string & qualities, vector<BamTools::CigarOp> & cigar_data) {

    double probability = 1;

    vector<unsigned short> deletions;
    vector<unsigned short> insertions;

    vector<BamTools::CigarOp>::iterator cigar_data_it = cigar_data.begin();

    unsigned short running_length = 0;

    while (cigar_data_it != cigar_data.end()) {

        if (cigar_data_it->Type == 'M') {

            running_length += cigar_data_it->Length;
        
        } else if (cigar_data_it->Type == 'D') {

            deletions.push_back(running_length);
        
        } else if (cigar_data_it->Type == 'I') {

            for (uint i = 0; i < cigar_data_it->Length; i++) {

                insertions.push_back(running_length);                
                running_length++;
            }
        }

        cigar_data_it++;
    }

    assert(running_length == qualities.size());

    deletions.push_back(qualities.size() + 1);
    insertions.push_back(qualities.size() + 1);

    vector<unsigned short>::iterator deletions_it = deletions.begin();
    vector<unsigned short>::iterator insertions_it = insertions.begin();

    boost::regex digit_expression("\\d+");
    boost::sregex_token_iterator iter(md_tag.begin(), md_tag.end(), digit_expression, 0);
    boost::sregex_token_iterator end;

    unsigned short q_idx = 0;
    unsigned short cigar_count = 0;

    while (iter != end) {

        for (unsigned short i = 0; i < atoi(iter->str().c_str()); i++) {

            while (*insertions_it == q_idx) {

                insertions_it++;
                q_idx++;
            }

            probability *= (1 - pow(10,-((int) qualities[q_idx] - 33)/(double)10));
            q_idx++;
        } 

        while (*insertions_it == q_idx) {

            insertions_it++;
            q_idx++;
        }

        iter++;

        if (q_idx < qualities.size()) {
            
            if (*deletions_it == q_idx) {
                
                deletions_it++;

            } else {

                probability *= pow(10,-((int) qualities[q_idx] - 33)/(double)10)/3;
                q_idx++;
            }
        
        } else {

            assert(iter == end);
        }
    }

    assert(*insertions_it == qualities.size() + 1);
    assert(*deletions_it == qualities.size() + 1);
    assert(q_idx == qualities.size());
    assert(probability < 1);
    assert(probability >= double_underflow);

    return probability;
}


void Assembler::bayesemblerCallback(deque< pair<list<string>, vector<FragmentAlignment *> > *> * graph_queue, deque<AssemblyInfo *> * output_queue, boost::mutex * graph_lock, boost::mutex * output_lock, int * bayesembler_remaining_graphs, double adjusted_fragment_count, FragmentLengthModel * fragment_length_model_pt) {

	while (true) {
 
        graph_lock->lock();

        if ((*bayesembler_remaining_graphs) == 0) {

            graph_lock->unlock();
            break;
        }

		if (graph_queue->size() > 0) {

			pair<list<string>, vector<FragmentAlignment *> > * current_graph_pt = graph_queue->front();
			graph_queue->pop_front();
			(*bayesembler_remaining_graphs)--;
            int graph_seed = options_variables.seed + (*bayesembler_remaining_graphs);

			graph_lock->unlock();

            AssemblyInfo * assembly = new AssemblyInfo;
            stringstream assembly_log_stream;
            
            boost::random::mt19937 mt_rng;
            mt_rng.seed(graph_seed);
            
            assembly_log_stream << "\n[" << getLocalTime() << "] " << "Seeding graph random number generator with " << graph_seed << endl;
                                    
            // Parse dot-file and tabulate all paths in the graph
            vector <Candidate> all_candidates;
            int max_transcript_length = 0;
            stringstream assembly_id_ss;
            int sub_graph_count = 0; 

            // Loop over dot-strings
            int candidate_counter = 0; 

            tr1::unordered_map<string, int> subgraph_complexity;
            
            for (list<string>::iterator dot_iter = current_graph_pt->first.begin(); dot_iter != current_graph_pt->first.end(); dot_iter++) {

                assembly_log_stream << "\n[" << getLocalTime() << "] " << "Calculating paths for subgraph " << sub_graph_count << endl;           

                AllPaths paths(*dot_iter);
                vector <Candidate> current_candidates = paths.findFullPaths(options_variables.no_pre_mrna, options_variables.max_candidate_number, assembly_log_stream);

                for (vector<Candidate>::iterator candidate_iter = current_candidates.begin(); candidate_iter != current_candidates.end(); candidate_iter++) {

                    // Reset of candidate index due to current alignmentparser implementation
                    Candidate temp_candidate = *candidate_iter; 
                    temp_candidate.idx = candidate_counter;
                    candidate_counter++;

                    all_candidates.push_back(temp_candidate);

                    if (temp_candidate.length > max_transcript_length) {
                    
                        max_transcript_length = temp_candidate.length;
                    }
                }
                
                assembly_id_ss << current_candidates.front().graph_name; 

                assert(subgraph_complexity.insert(pair<string, uint>(current_candidates.front().graph_name, current_candidates.size())).second);

                if (sub_graph_count < current_graph_pt->first.size() - 1) {

                    assembly_id_ss << " ";
                }

                sub_graph_count++;    
            }
            
            string assembly_name = assembly_id_ss.str();    
            assembly->name = assembly_name;

            if (all_candidates.empty()) {
                
                assembly->return_code = 2;
                assembly->assembly_log = assembly_log_stream.str();

                output_lock->lock();
                output_queue->push_back(assembly);
                output_lock->unlock();
                continue;   
            }

            assert(!subgraph_complexity.empty());

            // Init sequencing model and parse alignment
            SequencingModel * sequencing_model = new SequencingModel(fragment_length_model_pt, max_transcript_length, assembly_log_stream);
            
            AlignmentParser alignment_parser;            
            CollapsedMap frag_tran_map = alignment_parser.calculateFragTranProbabilities(all_candidates, &(current_graph_pt->second), sequencing_model, false, assembly_log_stream, &subgraph_complexity);
         
            delete sequencing_model;
            delete current_graph_pt;
            
            if (frag_tran_map.return_code != 0) {
                
                assembly->return_code = frag_tran_map.return_code;
                assembly->assembly_log = assembly_log_stream.str();

                // cout << "DEBUG: Thread " << boost::this_thread::get_id() << " does frag tran map exit with after alignmentparser on graph " << dot_string_prefix_debug.str() << endl;

    			output_lock->lock();
    			output_queue->push_back(assembly);
    			output_lock->unlock();
                continue;	
            }
                        
            int num_transcripts = frag_tran_map.probability_matrix.cols();
            assert (frag_tran_map.new_idx_to_old_idx_effective_length.size() == num_transcripts);
            int num_fragments = 0;
            
            for (int i=0; i < frag_tran_map.counts.size(); i++) {
                
                num_fragments += frag_tran_map.counts[i];
            }
            
            assert (num_transcripts > 0);
            assert (num_fragments > 0);

            assembly_log_stream << "\n[" << getLocalTime() << "] Number of candidates used for further analysis: " << num_transcripts << endl;
            assembly_log_stream << "[" << getLocalTime() << "] Number of fragments used for further analysis: " << num_fragments << endl;
            
            // Set assignment options
            AssignmentGenerator assignment_generator(&mt_rng, &frag_tran_map.probability_matrix, &frag_tran_map.counts);
                                   
            // Calculate the minimum transcript set needed to explain all reads                            
            assignment_generator.calcMinimumReadCover(assembly_log_stream);                
            double msc_to_all_ratio = double(assignment_generator.getMinimumReadCoverSize())/num_transcripts;

            if (msc_to_all_ratio > double_almost_one) {

                msc_to_all_ratio = double_almost_one;
            }

            assert(msc_to_all_ratio >= double_underflow);
                    
            // Set simplex size options
            assembly_log_stream << "[" << getLocalTime() << "] " << "Initialising Gibbs sampler simplex size generator with minimum read based fixed pi and fixed gamma (pi=" << msc_to_all_ratio << ", gamma=" << options_variables.gamma << ")" << endl;
            FixedBinomialFixedGammaSimplexSizeGenerator fixed_simplex_size_generator = FixedBinomialFixedGammaSimplexSizeGenerator(msc_to_all_ratio, options_variables.gamma, num_transcripts, &mt_rng, num_fragments);
            SimplexSizeGenerator * simplex_size_generator = &fixed_simplex_size_generator;
                                  
            // Set gamma sampling options                                
            assert (options_variables.gamma >= double_underflow);
            assembly_log_stream << "[" << getLocalTime() << "] " << "Initialising Gibbs sampler gamma generator with fixed gamma (gamma=" << options_variables.gamma << ")" << endl;
            FixedGammaGenerator fixed_gamma_generator = FixedGammaGenerator(options_variables.gamma);
            GammaGenerator * gamma_generator = &fixed_gamma_generator;
            
            // Set expression option
            assembly_log_stream << "[" << getLocalTime() << "] " << "Initialising Gibbs sampler expression generator with symmetric Dirichlet prior" << endl;
            ExpressionGenerator expression_generator(num_fragments, num_transcripts, &mt_rng);                                        
                             
            // Init Gibbs sampler
            GibbsSampler gibbs_sampler(simplex_size_generator, gamma_generator, &expression_generator, &assignment_generator, frag_tran_map.new_idx_to_old_idx_effective_length, adjusted_fragment_count, num_fragments, num_transcripts, &mt_rng, false);

            int max_subgraph_complexity = 0;

            for (tr1::unordered_map<string, int>::iterator it = subgraph_complexity.begin(); it != subgraph_complexity.end(); it++) {

                assert(it->second >= 0);
                max_subgraph_complexity = max(max_subgraph_complexity, it->second);
            }

            assert(max_subgraph_complexity > 0);

            if (subgraph_complexity.size() == 1) {

                assert(num_transcripts == max_subgraph_complexity);
            
            } else {

                assert(num_transcripts >= max_subgraph_complexity);
            }

            uint gibbs_burn_in = options_variables.gibbs_base_iterations + options_variables.gibbs_scale_iterations * max_subgraph_complexity;
            uint gibbs_samples = 10 * options_variables.gibbs_base_iterations + 10 * options_variables.gibbs_scale_iterations * max_subgraph_complexity;
            
            // Run sampler 
            GibbsOutput gibbs_output(gibbs_samples, num_transcripts);
            gibbs_sampler.runSampler(&gibbs_output, gibbs_burn_in, gibbs_samples, assembly_log_stream);

            assert (frag_tran_map.return_code == 0);
            assembly->return_code = 0;
                
            // Get marginal indices 
            Ensemble top_marginals = gibbs_output.getTopMarginals(options_variables.marginal_threshold);
            
            if (top_marginals.indices.size() < 1) {
            
                assembly_log_stream << "\n[" << getLocalTime() << "] No transcripts with a marginal probability over " << options_variables.marginal_threshold << " found - exiting!" << endl;                  
                assembly->return_code = 5;
                
            } else {
                        
                assembly->ensemble = gibbs_output.writeEnsembleGtf(top_marginals, all_candidates, frag_tran_map.new_idx_to_old_idx_effective_length, adjusted_fragment_count, options_variables.count_threshold);

                if (assembly->ensemble.size() < 1) {

                    assembly_log_stream << "\n[" << getLocalTime() << "] No transcripts with an expected count threshold over " << options_variables.count_threshold << " found - exiting!" << endl;                  
                    assembly->return_code = 7;
                }                       
            }    
            
            if (options_variables.output_mode == "full") {

                assembly->candidates = gibbs_output.writeCandidateGtf(all_candidates, frag_tran_map.new_idx_to_old_idx_effective_length, adjusted_fragment_count);
                assembly->assembly_log = assembly_log_stream.str();
            }
            
			// Process Bayesembler output
			output_lock->lock();
			output_queue->push_back(assembly);
			output_lock->unlock();

		} else {
		    
            graph_lock->unlock();
		}
	}
}


void Assembler::fragmentLengthCallback(deque< pair<list<string>, vector<FragmentAlignment *> > *> * graph_queue, map<uint,uint> * output_map, boost::mutex * graph_lock, boost::mutex * output_lock, int * remaining_graphs, uint * single_path_long_candidates) {

    map<uint,uint> fragment_length_map;
    uint cur_single_path_long_candidates = 0;

    while (true) {
        
        graph_lock->lock();

        if ((*remaining_graphs) == 0) {

            graph_lock->unlock();
            break;
        }

        if (graph_queue->size() > 0) {

            pair<list<string>, vector<FragmentAlignment *> > * current_graph_pt = graph_queue->front();
            graph_queue->pop_front();
            (*remaining_graphs)--;
            graph_lock->unlock();
                              
            stringstream fragment_length_log_stream;
                  
            // Parse dot-file and tabulate all paths in the graph
            vector <Candidate> all_candidates;
            stringstream assembly_id_ss;
            int sub_graph_count = 0; 

            // Loop over dot-strings
            int candidate_counter = 0; 

            assert (current_graph_pt->first.size() == 1);

            for (list<string>::iterator dot_iter = current_graph_pt->first.begin(); dot_iter != current_graph_pt->first.end(); dot_iter++) {

                fragment_length_log_stream << "\n[" << getLocalTime() << "] " << "Calculating paths for subgraph " << sub_graph_count << endl;           

                AllPaths paths(*dot_iter);
                vector <Candidate> current_candidates = paths.findFullPaths(true, options_variables.max_candidate_number, fragment_length_log_stream);
                
                for (vector<Candidate>::iterator candidate_iter = current_candidates.begin(); candidate_iter != current_candidates.end(); candidate_iter++) {

                    // Reset of candidate index due to current alignmentparser implementation
                    Candidate temp_candidate = *candidate_iter;

                    if (temp_candidate.length < options_variables.fragment_est_min_transcript_length) {

                        continue;
                    
                    } else {

                        cur_single_path_long_candidates++;
                    }

                    temp_candidate.idx = candidate_counter;
                    candidate_counter++;

                    all_candidates.push_back(temp_candidate);
                }

                sub_graph_count++;    
            }

            // Init sequencing model and parse alignment 

            if (!all_candidates.empty()) {         
            
                AlignmentParser alignment_parser;        
                alignment_parser.fetchFragmentLengths(&fragment_length_map, all_candidates, &(current_graph_pt->second), fragment_length_log_stream);
            }

            delete current_graph_pt;
                                
        } else {
            
            graph_lock->unlock(); 
        }
    }

    output_lock->lock();

    *single_path_long_candidates += cur_single_path_long_candidates;

    for (map<uint,uint>::iterator fragment_length_map_iter = fragment_length_map.begin(); fragment_length_map_iter != fragment_length_map.end(); fragment_length_map_iter++) {

        if (output_map->count(fragment_length_map_iter->first) > 0) {

            (*output_map)[fragment_length_map_iter->first] += fragment_length_map_iter->second; 

        } else {

            (*output_map)[fragment_length_map_iter->first] = fragment_length_map_iter->second; 
        }
    }

    output_lock->unlock();
}


void Assembler::outputCallback(deque<AssemblyInfo *> * output_queue, boost::mutex * output_lock, int * output_remaining_graphs, int num_graphs) {

    ofstream ensemble_file("assembly.gtf");
    ofstream candidates_file;
    ofstream log_file;


    if (options_variables.output_mode == "full") {

        candidates_file.open("candidates.gtf");
        log_file.open("graph_logs.txt");   
    }

    vector <int> stats_list(8,0);
    int bayesembler_counter = 0;
    
	while (true) {
        
        output_lock->lock();
        if ((*output_remaining_graphs) == 0) {

            output_lock->unlock();
            break;
        }

		if (output_queue->size() > 0) {
            // cout << "DEBUG: output_remaining_graphs " << *output_remaining_graphs << endl;	        
            AssemblyInfo * current_assembly_pt = output_queue->front();
			output_queue->pop_front();
			(*output_remaining_graphs)--;
            output_lock->unlock();
      
            if (current_assembly_pt->return_code == 0) {
            
                ensemble_file << current_assembly_pt->ensemble;
                stats_list[current_assembly_pt->return_code]++;

            } else if (current_assembly_pt->return_code == 5 or current_assembly_pt->return_code == 7) {
       
                stats_list[current_assembly_pt->return_code]++;                
            
            } else {
                
                assert (current_assembly_pt->candidates.size() == 0);
                stats_list[current_assembly_pt->return_code]++;  
            }
            
            bayesembler_counter++;

            if (options_variables.output_mode == "full") {
                
                candidates_file << current_assembly_pt->candidates;
                log_file << "\n###### " << current_assembly_pt->name << " ######\n" << current_assembly_pt->assembly_log << endl;
            }
            
            if (bayesembler_counter % 2000 == 0) {

                output_lock->lock();        
                cout << "[" << getLocalTime() << "] Completed assembly and post-processing of " << bayesembler_counter << " graph(s) with " << output_queue->size() << " in output queue and " << num_graphs-bayesembler_counter-(*output_remaining_graphs) << " skipped (" << *output_remaining_graphs << " remaining)" << endl;             
                output_lock->unlock();         
            }

            delete current_assembly_pt;
     
		} else {
		    
            output_lock->unlock();
            boost::this_thread::sleep(boost::posix_time::seconds(10));      
		}

	}
        
    cout << "\n[" << getLocalTime() << "] " << "======= FINAL STATS =======\n\n";
    cout << "\t" << num_graphs << " graph(s) were parsed in total" << endl;
    cout << "\t" << stats_list[0] << " graph(s) were assembled successfully" << endl;
    cout << "\t" << num_graphs-bayesembler_counter << " graph(s) were not assembled due to non successfully assigned fragments" << endl;
    cout << "\t" << stats_list[2] << " graph(s) were not assembled due to insufficient exon coverage" << endl;
    cout << "\t" << stats_list[3] << " graph(s) were not assembled due to insufficient candidate coverage" << endl;
    cout << "\t" << stats_list[4] << " graph(s) were not assembled as no fragments had sufficient mapping probability" << endl;
    cout << "\t" << stats_list[5] << " graph(s) were not assembled as no transcripts had probability higher than " << options_variables.marginal_threshold << endl;
    cout << "\t" << stats_list[7] << " graph(s) were not assembled due to transcripts being pre-mRNA or having an expected count lower than " << options_variables.count_threshold << endl;
    cout << endl;
    
    ensemble_file.close();

    if (options_variables.output_mode == "full") {

        candidates_file.close();
        log_file.close();
    }
}


void Assembler::runBayesemblerThreaded() {
    
	// Infer splice-graphs
    double adjusted_fragment_count = 0;
	vector<pair<list<GraphInfo>, string> > graph_list = generateSpliceGraphs(&adjusted_fragment_count);
    
    vector<pair<list<GraphInfo>, string> > single_graph_list; 
    vector<pair<list<GraphInfo>, string> > multi_graph_list;

    int single_graph_number = 0;
    int multi_graph_number = 0;

    cout << "[" << getLocalTime() << "] " << "Sorting splice-graphs by read count" << endl;
    
    // Separate out single-path graphs and sort multi-path graphs according to mapped read number
    for (int i=0; i < graph_list.size(); i++) {

        list<GraphInfo> single_graph_list_temp;
        list<GraphInfo> multi_graph_list_temp;
        list<pair<uint,uint> > multi_graph_list_idx_temp;
        uint graph_list_idx = 0;

        for (list<GraphInfo>::iterator graph_iter = graph_list[i].first.begin(); graph_iter != graph_list[i].first.end(); graph_iter++) {

            if (graph_iter->single_path_graph) {

                assert (graph_iter->dot_strings.size() == 1);
                single_graph_list_temp.push_back(*graph_iter);
                single_graph_number++;
            
            } else {

                assert (graph_iter->dot_strings.size() >= 1);

                multi_graph_list_temp.push_back(*graph_iter);
                multi_graph_number++;

                pair<uint,uint> temp_pair(graph_list_idx, graph_iter->read_count);
                multi_graph_list_idx_temp.push_back(temp_pair);
                graph_list_idx++;

            }
        }

        multi_graph_list_idx_temp.sort(Assembler::graphCompareReadCount);
        list<GraphInfo> multi_graph_list_temp_sorted;

        list<list<GraphInfo>::iterator> multi_graph_erase_iter_list;
        list<pair<uint,uint> >::iterator sort_idx_iterator = multi_graph_list_idx_temp.begin(); 

        // Fetch first 1000 graphs according to sorted index
        uint prio_sort_count = 0;
        // while (sort_idx_iterator != multi_graph_list_idx_temp.end() and prio_sort_count < 1000) {
        while (sort_idx_iterator != multi_graph_list_idx_temp.end()) {

            prio_sort_count++;

            list<GraphInfo>::iterator multi_graph_list_temp_iter = multi_graph_list_temp.begin();
            advance(multi_graph_list_temp_iter, sort_idx_iterator->first);

            multi_graph_list_temp_sorted.push_back(*multi_graph_list_temp_iter);
            multi_graph_erase_iter_list.push_back(multi_graph_list_temp_iter);
            sort_idx_iterator++;
        }

        assert (multi_graph_list_temp.size() >= multi_graph_erase_iter_list.size());

        for (list<list<GraphInfo>::iterator>::iterator graph_erase_iter = multi_graph_erase_iter_list.begin(); graph_erase_iter != multi_graph_erase_iter_list.end(); graph_erase_iter++) {
            
            multi_graph_list_temp.erase(*graph_erase_iter);

        }

        // Add the remaining graphs in chromosome order
        multi_graph_list_temp.splice(multi_graph_list_temp.begin(), multi_graph_list_temp_sorted);

        assert(multi_graph_list_temp_sorted.size() == 0);
        assert(graph_list[i].first.size() == single_graph_list_temp.size() + multi_graph_list_temp.size());


        pair<list<GraphInfo>, string> single_pair(single_graph_list_temp, graph_list[i].second);
        single_graph_list.push_back(single_pair);
 
        pair<list<GraphInfo>, string> multi_pair(multi_graph_list_temp, graph_list[i].second);
        multi_graph_list.push_back(multi_pair);

    }

    cout << "[" << getLocalTime() << "] " << "Finished sorting splice-graphs by read count\n" << endl;

    boost::mutex graph_lock;
    boost::mutex output_lock;

    // Estimate fragment lengths
    FragmentLengthModel * fragment_length_model_pt; 

    uint single_path_long_candidates = 0;
    
    if (options_variables.estimate_fragment_length) {

        vector<pair<list<GraphInfo>, string> > single_graph_list_fragment = single_graph_list;

        cout << "[" << getLocalTime() << "] Spawning " << options_variables.num_threads << " thread(s) for fetching alignments and 1 i/o thread" << endl;

        int single_remaining_graphs = single_graph_number;
        int single_output_remaining_graphs = single_graph_number;

        deque< pair<list<string>, vector<FragmentAlignment *> > *> fragment_length_graph_queue;
        map<uint,uint> output_map;

        boost::thread fragment_length_grapher_thread(boost::bind(&Assembler::grapherCallback, this, &single_graph_list_fragment, &fragment_length_graph_queue, &graph_lock, &output_lock, &single_remaining_graphs, &single_output_remaining_graphs, 10*options_variables.num_threads, true));
        boost::thread_group fragment_length_bayesembler_threads;
            
        for (int i = 0; i < options_variables.num_threads; i++) {
            
            fragment_length_bayesembler_threads.create_thread(boost::bind(&Assembler::fragmentLengthCallback, this, &fragment_length_graph_queue, &output_map, &graph_lock, &output_lock, &single_remaining_graphs, &single_path_long_candidates));
        }

        fragment_length_bayesembler_threads.join_all();
        fragment_length_grapher_thread.join();

        cout << "[" << getLocalTime() << "] Estimating fragment length distribution from " << single_path_long_candidates << " transcripts longer than " << options_variables.fragment_est_min_transcript_length << " nucleotides" << endl; 

        if (options_variables.output_mode == "full") {

            ofstream fragment_length_file("fragment_lengths.txt");
            fragment_length_file << "length\tcount" << endl;

            for (map<uint,uint>::iterator output_map_it = output_map.begin(); output_map_it != output_map.end(); output_map_it++) {

                fragment_length_file << output_map_it->first << "\t" << output_map_it->second << endl;
            }

            fragment_length_file.close();
        }

        fragment_length_model_pt = new GaussianFragmentLengthModel(output_map);
    
    } else {

        cout << "[" << getLocalTime() << "] Using user-defined fragment length distribution parameters." << endl;

        fragment_length_model_pt = new GaussianFragmentLengthModel(options_variables.frag_mean, options_variables.frag_sd);
    }

    int graph_number = single_graph_number + multi_graph_number; 

    cout << "\n[" << getLocalTime() << "] Starting Bayesembler on " << multi_graph_number << " multi-path graph(s) and " << single_graph_number << " single-path graph(s)" << endl;

    int bayesembler_remaining_graphs = graph_number;
    int output_remaining_graphs = graph_number; 

    vector<pair<list<GraphInfo>, string> > all_graphs_sorted;

    for (vector<pair<list<GraphInfo>, string> >::iterator multi_graph_list_iter = multi_graph_list.begin(); multi_graph_list_iter != multi_graph_list.end(); multi_graph_list_iter++) {

        all_graphs_sorted.push_back(*multi_graph_list_iter);
    }

    for (vector<pair<list<GraphInfo>, string> >::iterator single_graph_list_iter = single_graph_list.begin(); single_graph_list_iter != single_graph_list.end(); single_graph_list_iter++) {

        all_graphs_sorted.push_back(*single_graph_list_iter);
    }

    cout << "[" << getLocalTime() << "] Spawning " << options_variables.num_threads << " Bayesembler thread(s) and 2 i/o threads\n" << endl;

	deque< pair<list<string>, vector<FragmentAlignment *> > *> graph_queue;
    deque< AssemblyInfo *> output_queue;

    boost::thread grapher_thread(boost::bind(&Assembler::grapherCallback, this, &all_graphs_sorted, &graph_queue, &graph_lock, &output_lock, &bayesembler_remaining_graphs, &output_remaining_graphs, 10*options_variables.num_threads, options_variables.keep_temp_files));
	boost::thread output_thread(boost::bind(&Assembler::outputCallback, this, &output_queue, &output_lock, &output_remaining_graphs, graph_number));

	// Spawn Bayesembler worker threads
	boost::thread_group bayesembler_threads;
    	
	for (int i = 0; i < options_variables.num_threads; i++) {
    	
    	bayesembler_threads.create_thread(boost::bind(&Assembler::bayesemblerCallback, this, &graph_queue, &output_queue, &graph_lock, &output_lock, &bayesembler_remaining_graphs, adjusted_fragment_count, fragment_length_model_pt));
    }

	bayesembler_threads.join_all();
    grapher_thread.join();
	output_thread.join();

    delete fragment_length_model_pt;
}

