
/*
alignmentParser.cpp - This file is part of the Bayesembler (v1.1.1)


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


#include <utils.h>
#include <map>
#include <tr1/unordered_map>
#include <algorithm>
#include <set>
#include <alignmentParser.h>
#include <sequencingModel.h>
#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include <api/BamAux.h>
#include <boost/lexical_cast.hpp>

pair<bool, AlignmentParser::FragmentMatch> AlignmentParser::matchFragmentToCandidate(FragmentAlignment alignment, Candidate candidate, string strand) {

    vector<Vertex> current_vertices = candidate.vertices;

    int current_pos = alignment.read_data.first.first + 1;
    uint first_exon_idx = 0;
    uint current_exon_idx = 0;

    vector<BamTools::CigarOp> current_cigar = alignment.read_data.first.second;
    
    bool skip_fragment = true;

    // Loop over reads in pair
    for (uint pair_idx = 0; pair_idx < 2; pair_idx++) {

        bool first_match = true;

        if (pair_idx == 1) {

            current_pos = alignment.read_data.second.first + 1;
            current_cigar = alignment.read_data.second.second;
        }

        assert (current_cigar.front().Type == 'M');
        assert (current_cigar.back().Type == 'M');

        // Candidate matching
        for (uint j=0; j < current_cigar.size(); j++) {
            
            skip_fragment = true;
            
            // If match character
            if (current_cigar[j].Type == 'M') {
            
                // If first match in cigarstring, search for exon where alignment starts
                if (first_match) {
                                        
                    for (uint k = first_exon_idx; k < current_vertices.size(); k++) {
                        
                        if ((current_pos >= current_vertices[k].start) and (current_pos <= current_vertices[k].end)) {
                            
                            if (pair_idx == 0) {

                                first_exon_idx = k;
                            }
                                    
                            // Set position to end of match
                            current_pos += current_cigar[j].Length - 1;
                            current_exon_idx = k;
                            first_match = false;
                            skip_fragment = false;
                            break;
                        }
                    }    
                
                // If last character was not an intron, update current position and move on 
                } else if (current_cigar[j-1].Type != 'N') {
                                            
                    assert (current_cigar[j-1].Type != 'M');
                    
                    if (current_pos <= current_vertices[current_exon_idx].end) {

                        current_pos += current_cigar[j].Length - 1;
                        skip_fragment = false;  
                    }
                
                // If not first match, check that current position matches the next exon start (splice acceptor site check)
                } else {
                    
                    assert (current_cigar[j-1].Type != 'M');
                    
                    if (current_vertices[current_exon_idx].start == current_pos) {
                        
                        current_pos += current_cigar[j].Length - 1;
                        skip_fragment = false;
                    }                 
                }
            
            // If ref-skip character, check donor exon-intron boundary (splice donor site check), increment exon idx and update current position to next exon start
            } else if (current_cigar[j].Type == 'N') {
            
                if (current_vertices[current_exon_idx].end == current_pos) {
                    
                    current_exon_idx += 1;
                    current_pos += current_cigar[j].Length + 1;

                    if (current_exon_idx >= current_vertices.size()) {

                        skip_fragment = true;
                    
                    } else {

                        skip_fragment = false;
                    }
                }
            
            // IMPORTANT: The following handling of insertions and deletions contradict the samformat specification, but conforms with TopHat's definitions.
            // If insertion to the reference, no action
            } else if (current_cigar[j].Type == 'I') {
                
                if (current_pos <= current_vertices[current_exon_idx].end) {
                
                    skip_fragment = false;            
                }
            
            // If deletion in the reference, increment position
            } else if (current_cigar[j].Type == 'D') {
                
                if (current_pos <= current_vertices[current_exon_idx].end) {

                    current_pos += current_cigar[j].Length + 1;
                    skip_fragment = false;
                }
                
            } else {
                
                cerr << "ERROR: Unhandled cigar string symbol '" << current_cigar[j].Type << "'!" << endl;
                exit(-1);
            }
                
            if (skip_fragment) {
        
                break;
            }
            
            assert(current_exon_idx >= 0);
        }
     
        // Final right boundary check
        if (skip_fragment == false) {
                            
            assert (first_exon_idx >= 0);
            assert (current_exon_idx >= first_exon_idx);
            assert(current_exon_idx < current_vertices.size());
            
            if (current_pos > current_vertices[current_exon_idx].end) {
            
                skip_fragment = true;                
            }
        } 

        if (skip_fragment) {
    
            break;
        }
    }

    FragmentMatch fragment_match; 

    if (!skip_fragment) {

        fragment_match.first_exon_idx = first_exon_idx;
        fragment_match.last_exon_idx = current_exon_idx;

        // Calculate fragment length
        if (first_exon_idx == current_exon_idx) {

            fragment_match.fragment_length = current_pos - (alignment.read_data.first.first + 1) + 1;

        } else {

            fragment_match.fragment_length = current_vertices[first_exon_idx].end - (alignment.read_data.first.first + 1) + 1;

            for (uint j= (first_exon_idx +1); j < current_exon_idx; j++) {

                fragment_match.fragment_length += current_vertices[j].end - current_vertices[j].start + 1;
            }

            fragment_match.fragment_length += current_pos - current_vertices[current_exon_idx].start + 1;
        }

        assert(fragment_match.fragment_length <= candidate.length);

        // If plus strand, take last transcript position of second read
        if (strand == "+") {

            fragment_match.three_prime_position = current_pos - current_vertices[current_exon_idx].start + 1;
            assert(fragment_match.three_prime_position >= 0);

            for (uint j=0; j < current_exon_idx; j++) {

                fragment_match.three_prime_position += current_vertices[j].end - current_vertices[j].start + 1;
            }
                                
        // If minus strand, take first transcript position of first read    
        } else {

            assert(strand == "-");
            fragment_match.three_prime_position = current_vertices[first_exon_idx].end - (alignment.read_data.first.first + 1) + 1;

            assert(fragment_match.three_prime_position >= 0);

            for (uint j=(first_exon_idx + 1); j < current_vertices.size(); j++) {

                fragment_match.three_prime_position += current_vertices[j].end - current_vertices[j].start + 1;
            }                                
        }

        assert(fragment_match.three_prime_position <= candidate.length);
        assert(fragment_match.three_prime_position >= 0);

        return pair<bool,FragmentMatch> (true, fragment_match);
    
    } else {

        return pair<bool,FragmentMatch> (false, fragment_match);        
    }
}

CollapsedMap AlignmentParser::calculateFragTranProbabilities(vector<Candidate> &candidates, vector<FragmentAlignment*>* fragment_alignments, SequencingModel * sequencing_model, bool no_coverage_exclude, stringstream& log_stream, tr1::unordered_map<string, int> * subgraph_complexity) {
    
    uint row_buffer_size = 500;
    uint row_idx = 0;
    
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> probability_matrix(row_buffer_size, candidates.size());
    vector<uint> collapsed_counts;
        
    tr1::unordered_set <int> exclude_set;
	tr1::unordered_map <int, int> exclude_downstream_exons;
    tr1::unordered_map <int, int> exclude_upstream_exons;

    vector<uint> candidate_counts(candidates.size(), 0);
    
    for (uint i=0; i < candidates.size(); i++) {
                
        exclude_upstream_exons[i] = 2e9;
        exclude_downstream_exons[i] = -1;    
    }
        
    string first_strand = candidates.front().vertices.front().strand;
        	
    log_stream << "[" << getLocalTime() << "] " << "Matching fragments with candidates" << endl;
    
    // Stats
    uint matched_fragments = 0;
    uint unmatched_fragments = 0;
    uint fragment_candidate_matches = 0;
    uint underflow_skipped = 0;
    uint internal_exclude = 0;
    uint no_map_exclude = 0;
    uint upstream_exclude = 0;
    uint downstream_exclude = 0;
    
    for (vector<FragmentAlignment*>::iterator aligment_iter = fragment_alignments->begin(); aligment_iter != fragment_alignments->end(); aligment_iter++) {
    
        FragmentAlignment current_alignment = *(*aligment_iter); 

        bool found_fragment_match = false;
        
        // Init fragment probability vector
        Eigen::RowVectorXd temp_fragment_probs = Eigen::VectorXd::Zero(candidates.size());
        double norm_const = 0;
                         
        for (uint i=0; i < candidates.size(); i++) {
            
            Candidate current_candidate = candidates[i];

            if (exclude_set.count(i)) {
                
                continue;
            }
                        
            // Match fragment to candidate
            pair<bool, FragmentMatch> fragment_match = matchFragmentToCandidate(current_alignment, current_candidate, first_strand);

            // If match found
            if (fragment_match.first) {

                found_fragment_match = true;
                fragment_candidate_matches++;  
                candidate_counts[i]++;
 
                double fragment_probability = 1;

                fragment_probability *= current_alignment.probability;
                fragment_probability *= sequencing_model->calculateThreePrimeProb(fragment_match.second.three_prime_position, current_candidate.length);
                fragment_probability *= sequencing_model->calculateLengthProb(fragment_match.second.fragment_length, fragment_match.second.three_prime_position);

                // Add fragment probability to temporary container
                temp_fragment_probs(i) = fragment_probability;
                norm_const += fragment_probability;
                                        
                // Check fragment coverage exclusion criterion (ignoring pre-mRNA)
                if (!no_coverage_exclude) {                  
    
                    // If fragment overlap, update downstream flank
                    if ((fragment_match.second.first_exon_idx <= exclude_downstream_exons[i]) or (candidate_counts[i] == 1)) {
                        
                        if (fragment_match.second.last_exon_idx > exclude_downstream_exons[i]) {

                            exclude_downstream_exons[i] = fragment_match.second.last_exon_idx;
                        }                           
        
                        if (fragment_match.second.first_exon_idx < exclude_upstream_exons[i]) {
            
                            exclude_upstream_exons[i] = fragment_match.second.first_exon_idx;
                        }

                    // If no overlap, exclude candidate
                    } else if (!current_candidate.isPreMrna) {
            
                        assert(exclude_downstream_exons[i] > -1);
                        assert(exclude_upstream_exons[i] < 2e9);
        
                        exclude_set.insert(i);
                        internal_exclude += 1;
    
                    } else {
        
                        assert(false);                  
                    }
                }
            }
        } // End candidate loop
        
        if (found_fragment_match) {
            
            matched_fragments++;

            if (norm_const < double_underflow) {
            
                underflow_skipped++;
                delete *aligment_iter;
                continue;
            }
                                                
            // Normalise and search for row matches
            temp_fragment_probs = temp_fragment_probs/norm_const;
            
            uint i;
            bool row_match = false;
            for (i=0; i < row_idx; i++) {
    
                row_match = true;
                for (uint j=0; j < probability_matrix.cols(); j++) {
        
                    if (!(double_compare(probability_matrix(i,j), temp_fragment_probs(j)))) {
            
                        row_match = false;
                        break;
                    }    
                }
    
                if (row_match) {
        
                    break;
                }  
            }

            if (row_match) {
    
                assert(i < collapsed_counts.size());
                collapsed_counts[i]++;
                pair<int,double> temp_pair(i, norm_const);
                                
            } else {
                
                if (row_idx == probability_matrix.rows()) {

                    Eigen::NoChange_t no_change;
                    probability_matrix.conservativeResize(probability_matrix.rows() + row_buffer_size, no_change);   
                }
                
                probability_matrix.row(row_idx) = temp_fragment_probs;
                collapsed_counts.push_back(1);
                pair<int,double> temp_pair(row_idx, norm_const);
                                
                row_idx++;
                assert(row_idx == collapsed_counts.size());
            }                 
        
        } else {
                    
            unmatched_fragments++;
        }

        delete *aligment_iter;

    } // End fragment loop
    
    assert(fragment_alignments->size() == matched_fragments + unmatched_fragments);

    // Exclusion on flank coverage and total fragment count (ignoring pre-mRNA)
    if (!no_coverage_exclude) {
        
        for (int i=0; i < candidates.size(); i++) {
                
            if (exclude_set.count(i) == 0) {
            
                // If no fragments mapped, exclude
                if (exclude_downstream_exons[i] == -1) {
                
                    no_map_exclude++;

                    // cout << "DEBUG: i " << i << " up " << exclude_upstream_exons[i] << " down " << exclude_downstream_exons[i] << " candidate counts " << candidate_counts << endl;
                    assert(candidate_counts[i] == 0);
                    assert(exclude_upstream_exons[i] == 2e9);
                    exclude_set.insert(i);

                    // cout  << "DEBUG: No fragment mapped to candidate " << current_idx << endl;

                    if (row_idx > 0) {
                        
                        assert (probability_matrix.block(0, i, row_idx, 1).sum() < double_underflow);
                    }
                            
                } else if (!candidates[i].isPreMrna) {

                    assert(exclude_upstream_exons[i] < 2e9);
                    assert(exclude_downstream_exons[i] > -1);
                    assert(exclude_upstream_exons[i] >= 0);
                    assert(exclude_downstream_exons[i] <= (candidates[i].vertices.size()-1));
                                
                    if (exclude_upstream_exons[i] > 0) {
                    
                        upstream_exclude += 1;
                        exclude_set.insert(i);
                                   
                    } else if (exclude_downstream_exons[i] < (candidates[i].vertices.size()-1)) {

                        downstream_exclude += 1;
                        exclude_set.insert(i);
                    }
                }
            }
        }

    } else {
        
        assert (exclude_set.size() == 0);
    }
    
    CollapsedMap collapsed_map;
          
    if (exclude_set.size() == candidates.size()) {
    
        log_stream << "\n[" << getLocalTime() << "] No candidates had sufficient fragment coverage to continue - exiting!" << endl;
        
        collapsed_map.return_code = 3;
        return collapsed_map;
    }

    log_stream << "[" << getLocalTime() << "] Excluding " << exclude_set.size() << " of " << candidates.size() << " candidates and renormalising conditional probability matrix " << endl;
    
    // Remove exluded candidates from matrix and generate idx map
    tr1::unordered_map<int, pair<int, double> > new_idx_to_old_idx_effective_length;
    int new_idx = 0;  

    vector <int>::iterator exclude_set_iter;    
    vector <int> exclude_vec (exclude_set.size(),0);
    
    int current_exclude_set_idx;
    
    if (exclude_set.size() > 0) {
    
        int exclude_vec_counter = 0;

        for (tr1::unordered_set<int>::iterator it = exclude_set.begin(); it != exclude_set.end(); it++) {
            
            assert(subgraph_complexity->count(candidates.at(*it).graph_name) > 0);
            (*subgraph_complexity)[candidates.at(*it).graph_name]--;

            exclude_vec[exclude_vec_counter] = *it;
            exclude_vec_counter++;
        }
    
        sort (exclude_vec.begin(), exclude_vec.end());
    
        exclude_set_iter = exclude_vec.begin();
        current_exclude_set_idx = *exclude_set_iter;
            
        assert (exclude_vec.size() == exclude_set.size());

    } else {
        
        current_exclude_set_idx = -1;   
    }
    
    for (int i=0; i < candidates.size(); i++) {
        
        assert(candidates[i].idx == i);
        
        if (i == current_exclude_set_idx) {
            
            // cout << "DEBUG: exclude idx " << i << endl;
            
            exclude_set_iter++;
            if (exclude_set_iter == exclude_vec.end()) {
                
                int block_cols = probability_matrix.cols() - current_exclude_set_idx - 1;
                
                if (block_cols > 0) {
                    
                    probability_matrix.block(0, new_idx, row_idx, block_cols) = probability_matrix.block(0, current_exclude_set_idx+1, row_idx, block_cols);
                }   

                current_exclude_set_idx = -1; 
                
            } else {
                
                int block_cols = (*exclude_set_iter) - current_exclude_set_idx - 1;
                
                assert((*exclude_set_iter) > current_exclude_set_idx);
                
                if (block_cols > 0) {
                                    
                    probability_matrix.block(0, new_idx, row_idx, block_cols) = probability_matrix.block(0, current_exclude_set_idx+1, row_idx, block_cols);
                }
                
                current_exclude_set_idx = *exclude_set_iter; 
            }
            
        } else {    
        
            pair<int,double> map_pair(i, sequencing_model->getEffectiveLength(candidates[i].length));    
            new_idx_to_old_idx_effective_length[new_idx] = map_pair;
            new_idx++;
        }            
    }
        
    assert(exclude_set.size() == no_map_exclude + internal_exclude + upstream_exclude + downstream_exclude);
    assert(new_idx == probability_matrix.cols()-exclude_set.size());
    assert(new_idx_to_old_idx_effective_length.size() == probability_matrix.cols()-exclude_set.size());
    
    // cout << "DEBUG: Probability matrix before final normalisation \n" << probability_matrix << endl;
    
    // Final normalisation (to correct for excluded candidates)            
    int num_cols = probability_matrix.cols() - exclude_set.size();
    vector<int> row_exclude;
        
    if (exclude_set.size() > 0) {
    
        for (int i=0; i < row_idx; i++) {
        
            double row_sum = probability_matrix.block(i, 0, 1, num_cols).sum();          
        
            // Remove row if all probabilities are smaller than almost zero 
            if (row_sum < double_underflow) {
                            
                underflow_skipped += collapsed_counts[i];
                row_exclude.push_back(i);
                        
            } else {

                Eigen::RowVectorXd norm_row = probability_matrix.block(i, 0, 1, num_cols)/row_sum;
            
                int j;
                bool row_match = false;
                for (j=0; j < i; j++) {

                    row_match = true;
                    for (int k=0; k < num_cols; k++) {
    
                        if (!(double_compare(probability_matrix(j,k), norm_row(k)))) {
        
                            row_match = false;
                            break;
                        }    
                    }

                    if (row_match) {
    
                        break;
                    }  
                }

                if (row_match) {
            
                    collapsed_counts[j] += collapsed_counts[i];
                    row_exclude.push_back(i);

                } else {
                            
                    probability_matrix.block(i, 0, 1, num_cols) = norm_row;
                }        
            }
        }
    }
    
    // Final output container 
    int num_rows = row_idx - row_exclude.size();
    collapsed_map.probability_matrix = Eigen::MatrixXd(num_rows, num_cols);
    collapsed_map.counts = vector<int>(num_rows);
    int final_row_idx = 0; 
    int exclude_idx = 0;

    for (int i=0; i < row_idx; i++) {
        
        if ((row_exclude.size() > 0) and (exclude_idx < row_exclude.size())) { 
        
            if (row_exclude[exclude_idx] != i) {
                
                collapsed_map.probability_matrix.row(final_row_idx) = probability_matrix.block(i, 0, 1, num_cols);
                collapsed_map.counts[final_row_idx] = collapsed_counts[i];
                final_row_idx++;

            } else {
                
                exclude_idx++;
            }
            
        } else {
            
            collapsed_map.probability_matrix.row(final_row_idx) = probability_matrix.block(i, 0, 1, num_cols);
            collapsed_map.counts[final_row_idx] = collapsed_counts[i];
            final_row_idx++;            
        }    
    }
    
    collapsed_map.new_idx_to_old_idx_effective_length = new_idx_to_old_idx_effective_length;
    
    assert (collapsed_map.probability_matrix.rows() == collapsed_map.counts.size());
    assert (collapsed_map.probability_matrix.rows() == num_rows);
    assert (final_row_idx == num_rows);
                
    if (collapsed_map.probability_matrix.rows() == 0) {
    
        log_stream << "\n[" << getLocalTime() << "] No rows in conditional probability matrix - exiting!" << endl;
        
        collapsed_map.return_code = 4;
        return collapsed_map;   
    }

    assert(exclude_set.size() == no_map_exclude + internal_exclude + upstream_exclude + downstream_exclude);
       
    // Print stats
    log_stream << "\n[" << getLocalTime() << "] === BAM-FILE PARSING STATS ===\n" << endl;
    log_stream << "  " << fragment_alignments->size() << " fragment(s) were parsed in total." << endl;
        
    log_stream << "\n  - Fragment-candidate match stats:" << endl;
    log_stream << "  " << matched_fragments << " fragment(s) matched a transcript." << endl;     
    log_stream << "  " << unmatched_fragments << " fragment(s) dit not match a transcript." << endl;     
    log_stream << "  " << underflow_skipped << " fragments were removed due to probability underflow." << endl;
    log_stream << "  " << fragment_candidate_matches << " fragment-candidate match(es) were found." << endl;     
    log_stream << "  " << collapsed_map.probability_matrix.rows() << " unique rows in probability matrix." << endl;
    
    log_stream << "\n  - Candidate exclusion stats:" << endl;
    log_stream << "  " << exclude_set.size() << " candidates were excluded due to low fragment coverage." << endl;
    log_stream << "  " << no_map_exclude << " candidates were excluded due to lack of mapped fragments." << endl;
    log_stream << "  " << internal_exclude << " candidates were excluded due to lack of internal fragment coverage." << endl;
    log_stream << "  " << upstream_exclude << " candidates were excluded due to lack of upstream fragment coverage." << endl;
    log_stream << "  " << downstream_exclude << " candidates were excluded due to lack of downstream fragment coverage." << endl;
    log_stream << endl;

    collapsed_map.return_code = 0;
    
    return collapsed_map;
}


void AlignmentParser::fetchFragmentLengths(map<uint,uint> * fragment_length_map, vector<Candidate> & candidates, vector<FragmentAlignment*>* fragment_alignments, stringstream& log_stream) {
            
    string first_strand = candidates.front().vertices.front().strand;

    log_stream << "[" << getLocalTime() << "] " << "Matching fragments with candidates" << endl;
        
    for (vector<FragmentAlignment*>::iterator aligment_iter = fragment_alignments->begin(); aligment_iter != fragment_alignments->end(); aligment_iter++) {
    
        FragmentAlignment current_alignment = *(*aligment_iter); 

        for (uint i=0; i < candidates.size(); i++) {
                
            Candidate current_candidate = candidates[i];
        
            // Match fragment to candidate
            pair<bool, FragmentMatch> fragment_match = matchFragmentToCandidate(current_alignment, current_candidate, first_strand);

            // If match found
            if (fragment_match.first) {
 
                if (fragment_length_map->count(fragment_match.second.fragment_length)) {

                    (*fragment_length_map)[fragment_match.second.fragment_length]++;

                } else {

                    (*fragment_length_map)[fragment_match.second.fragment_length] = 1;
                }
            }
        }

        delete *aligment_iter;
    }
}

