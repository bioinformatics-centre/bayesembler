
/*
main.cpp - This file is part of the Bayesembler (v1.2.0)


The MIT License (MIT)

Copyright (c) 2015 Lasse Maretty and Jonas Andreas Sibbesen

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


#include <string>
#include <time.h>
#include <boost/random.hpp>
#include <utils.h>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp> 
#include <assembler.h>
#include <unistd.h>

namespace po = boost::program_options;

typedef boost::random::mt19937* mt_rng_pt_t;

using namespace std;
				
int main (int argc, char * const argv[]) {
    
    OptionsContainer options_variables;
    
    cout << "\nYou are using the Bayesembler v1.2.0. For more information go to bayesembler.binf.ku.dk\n" << endl;

    try {
        po::options_description general("", 150);
        general.add_options()
            ("help,h", "produce help message for base options")
            ("advanced,a", "produce help message for advanced options")
		;
        
        po::options_description required("== Required arguments ==", 150);
        required.add_options()
            ("bam-file,b", po::value<string>(&options_variables.bam_file)->required(), "TopHat2 bamfile containing the mapped paired-end reads to the reference genome")      
        ;

		po::options_description optional("== Optional arguments ==", 150);
        optional.add_options()
            ("output-prefix,o", po::value<string>(&options_variables.output_prefix)->default_value(""), "prefix to be used on all output filenames (OBS: to use output directory <dir> without any prefix on the output files, use <dir/>)")
            ("num-threads,p", po::value<int>(&options_variables.num_threads)->default_value(1), "number of threads used for assembly (actual number of threads: <num-threads> + two I/O threads) ")
            ("strand-specific,s", po::value<string>(&options_variables.strand_specific)->default_value("unstranded"), "data is strand-specific. Use \"first\" to indicate mate orientation as in the dUTP protocol or \"second\" if opposite")
            ("confidence-threshold,c", po::value<double>(&options_variables.marginal_threshold)->default_value(0.5, "0.5"), "exclude candidates with a confidence below <confidence-threshold>")                        
            ("count-threshold,f", po::value<double>(&options_variables.count_threshold)->default_value(12,"12"), "exclude candidates with an expected fragment count below <count-threshold>")
            ("no-pre-mRNA,m", po::value<bool>(&options_variables.no_pre_mrna)->default_value(false)->implicit_value(true), "do not include pre-mRNA in the set of candidates")
        ;

        po::options_description advanced("== Advanced arguments ==", 150);
        advanced.add_options()
            ("seed", po::value<int>(&options_variables.seed)->default_value(time(NULL), "unix time"), "seed for pseudo-random number generator (produces reproducible results only when running with one thread)")
            ("output-mode", po::value<string>(&options_variables.output_mode)->default_value("assembly"), "output a gtf file of the assembly. Use \"full\" to output additional information on the assembly, candidates and fragment lengths")
            ("library-size", po::value<int>(&options_variables.total_num_fragments)->default_value(0), "the total number of sequenced paired-end reads used for FPKM normalisation. Use value of 0 (default) to normalise using the number of mapped paired-end reads")
            ("keep-temp-files", po::value<bool>(&options_variables.keep_temp_files)->default_value(false)->implicit_value(true), "keep filtered bam, instance and processsam log-file")
            ("max-candidate-number", po::value<int>(&options_variables.max_candidate_number)->default_value(100), "maximum number of candidates per splice-graph (used for graph pruning)")
            ("dirichlet-parameter", po::value<double>(&options_variables.gamma)->default_value(1), "abundance parameter (symmetric Dirichlet concentration parameter gamma)")
            ("frag-mean", po::value<double>(&options_variables.frag_mean), "disable internal fragment length mean estimation and use <frag-mean>")
            ("frag-sd", po::value<double>(&options_variables.frag_sd), "disable internal fragment length sd estimation and use <frag-sd>")
            ("gibbs-iteration-scaling", po::value<unsigned int>(&options_variables.gibbs_scale_iterations)->default_value(60), "scaling factor used to calculate number of Gibbs iterations (burn-in & sample-size) using: burn-in = 1000 + <gibbs-iteration-scaling> * number_of_candidates, sample-size = 10 * burn-in")    
		;

        /*	
			== PRESET PARAMETERS ==
        */

        // Lower bound on the length of transcripts used for fragment length estimation
        options_variables.fragment_est_min_transcript_length = 2500;

        // Number of Gibbs iterations (burn-in & sample-size) is calculated using: burn-in = <gibbs-base> + <gibbs-scale> * number_of_candidates, sample-size = 10 * burn-in"
		options_variables.gibbs_base_iterations = 1000;
    
		// Minimum fraction of bases with a non-zero coverage. Necessary for filtering out a few very long exons with no coverage assembled by the CEM splice-graph engine
        options_variables.exon_base_coverage = 0.05;

        po::options_description desc("## Bayesembler help ##");
        desc.add(general).add(required).add(optional);

        po::options_description adv("## Bayesembler help ##");
        adv.add(general).add(required).add(optional).add(advanced);
                
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, adv), vm);
        
        if (vm.count("help") || argc == 1) {
            cout << desc << endl;
            return 1;
        }
        
        if (vm.count("advanced")) {
            cout << adv << endl;
            return 1;
        }

        po::notify(vm);

        if (options_variables.output_prefix.size() > 0) {

            boost::filesystem::path prefix(options_variables.output_prefix);
            boost::filesystem::path prefix_dir = prefix.parent_path();

            if (!prefix_dir.empty()) {
                
                cout << "[" << getLocalTime() << "] Creating output directory " << prefix_dir << endl;
                boost::filesystem::create_directories(prefix_dir);
            }
        }

        if (options_variables.strand_specific != "unstranded" and options_variables.strand_specific != "first" and options_variables.strand_specific != "second") {
                
            cout << "\nERROR: --strand-specific argument need to be either \"unstranded\", \"first\" or \"second\"" << endl;
            exit(-1);
        }

        if (options_variables.output_mode != "assembly" and options_variables.output_mode != "full") {
                
            cout << "\nERROR: --output_mode argument need to be either \"assembly\" or \"full\"" << endl;
            exit(-1);
        }

        if (vm.count("frag-mean") or vm.count("frag-sd")) {

            assert(vm.count("frag-mean") and vm.count("frag-sd"));

            options_variables.estimate_fragment_length = false;
            cout << "[" << getLocalTime() << "] Using fragment length mean " << options_variables.frag_mean << " and sd " << options_variables.frag_sd << endl; 

        } else {

            options_variables.estimate_fragment_length = true;
        }

        // Getting path to executable this way only works on Linux
        char exe_path_buff[1024];
        ssize_t exe_path_len = readlink("/proc/self/exe", exe_path_buff, sizeof(exe_path_buff)-1);
        
        if (exe_path_len != -1) {
            
            exe_path_buff[exe_path_len] = '\0';
            string test = string(exe_path_buff);
        } 

        char * samtools_path_env;
        samtools_path_env = getenv("SAMTOOLS_PATH");
        if (samtools_path_env != NULL) {

            options_variables.samtools_path = string(samtools_path_env);
        
        } else if (exe_path_len != -1) {

            boost::filesystem::path exe_path(exe_path_buff);
            boost::filesystem::path exe_dir = exe_path.parent_path();
            boost::filesystem::path samtools_path = exe_dir / boost::filesystem::path("dependencies/samtools_0.1.19/samtools"); 
            options_variables.samtools_path = samtools_path.string();
        } 
        
        char * cem_processsam_path_env;
        cem_processsam_path_env = getenv("CEM_PROCESSSAM_PATH");
        if (cem_processsam_path_env != NULL) {

            options_variables.cem_processsam_path = string(cem_processsam_path_env);
        
        } else if (exe_path_len != -1) {

            boost::filesystem::path exe_path(exe_path_buff);
            boost::filesystem::path exe_dir = exe_path.parent_path();
            boost::filesystem::path cem_processsam_path = exe_dir / boost::filesystem::path("dependencies/cem/processsam"); 
            options_variables.cem_processsam_path = cem_processsam_path.string();
        }

        Assembler assembler(options_variables);
        assembler.runBayesemblerThreaded();     
	}
    
    catch(std::exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        return 1;
    }
    
    catch(...) {
        cerr << "Exception of unknown type!" << endl;
        return 1;
    }
    
    return 0;    
}
