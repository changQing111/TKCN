#include "kmc_file.h"
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <set>
std::string rc(std::string s) {
    std::reverse(s.begin(),s.end());
    for(char &c:s) c=c=='A'?'T':c=='T'?'A':c=='C'?'G':'C';
    return s;
}
int main(int argc,char**argv) {
 try {
    if(argc!=3)throw std::runtime_error("Usage: query_mmap DB_PREFIX markers.fa");
    CKMCFile db;
    if(!db.OpenForRA(argv[1]))throw std::runtime_error("Cannot open database");
    CKMCFileInfo info;
    if(!db.Info(info) || info.kmer_length!=31 || !info.both_strands || info.mode!=0)
        throw std::runtime_error("Requires canonical 31-mer occurrence database");
    std::ifstream in(argv[2]);if(!in)throw std::runtime_error("Cannot open query FASTA");
    std::string line,id,seq;std::set<std::string> seen;
    std::cout<<"kmer_id\tsequence\tstored_count\tstored_min_count\tcounter_bytes\n";
    auto emit=[&]() {
        if(id.empty())return;
        if(seq.size()!=31 || seq.find_first_not_of("ACGT")!=std::string::npos || !seen.insert(id).second)
            throw std::runtime_error("Invalid marker");
        std::string canonical=std::min(seq,rc(seq));CKmerAPI k(31);
        if(!k.from_string(canonical))throw std::runtime_error("Cannot encode kmer");
        uint64 count=0;if(!db.CheckKmer(k,count))count=0;
        std::cout<<id<<'\t'<<seq<<'\t'<<count<<'\t'<<info.min_count<<'\t'<<info.counter_size<<'\n';
    };
    while(std::getline(in,line)) {
        if(!line.empty() && line.back()=='\r')line.pop_back();
        if(line.empty())continue;
        if(line[0]=='>'){emit();id=line.substr(1);id=id.substr(0,id.find_first_of(" \t"));seq.clear();}
        else seq+=line;
    }
    emit();if(seen.empty())throw std::runtime_error("No queries");db.Close();return 0;
 } catch(const std::exception&e) {std::cerr<<e.what()<<'\n';return 1;}
}
