#include <iostream>
#include <stdio.h>
#include <stdlib.h> 
#include <vector>
#include <fstream>
#include <ctime>
#include <sstream>
#include <map>
#include <unordered_map>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iomanip>
#include "word_embedding.h"

//#include <boost/unordered_map.hpp> // unorder_map

using namespace std;

typedef unordered_map<string,int>::const_iterator testit;
typedef unordered_map<int,double>::const_iterator map_iter;


int main(int argc, char** argv) {
	
	Embedding embedding;
	embedding.display_attributes();
	embedding.load_Word2VecBinFormat("data/embeddings/GoogleNews-vectors-negative300.bin");
	embedding.display_attributes();
	float* temp;
	clock_t begin = clock();
	
	for(unsigned int i = 0 ; i < embedding.size_voc() ; i++){

		temp = embedding.get(embedding[i]);
		if(temp == nullptr){cout<<"Erreur mot non existant : "<< embedding[i] <<endl;}

	}


	clock_t end = clock();
  	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<< "Time to check the voc : "<< elapsed_secs <<endl;


	cout<<"Cosine similarity between test and tests : "<<embedding.cosine("test","tests") <<endl;

	cout<<"Cosine similarity between test and tests : "<<embedding.cosine("test","dfghtgh") <<endl;

	if(embedding.get("test") != nullptr ){

		cout<<embedding.get("test")<<endl;

	}

	if(embedding.get("dfghtgh") != nullptr ){

		cout<<embedding.get("dfghtgh")<<endl;

	}
	
	

/*
	begin = clock();

	results = basic_language_model(queries , collection , k);
	
	end = clock();
  	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<< "Time to perform all the queries : "<< elapsed_secs <<endl;

	write_res_file(results , "data/res/results_basic" , "CHIC-");

	Embedding embedding;
	embedding.display_attributes();

	cout<<"loading "<< file_name.c_str() <<endl;

	clock_t begin = clock();

	embedding.load_Word2VecBinFormat(file_name.c_str());

	clock_t end = clock();
  	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<< "Time to read the embedded vectors : "<< elapsed_secs <<endl;

	embedding.display_attributes();

	cout<<"Cosine similarity between test and tests : "<<embedding.cosine("test","tests") <<endl;

	cout<<"Cosine similarity between test and tests : "<<embedding.cosine("test","dfghtgh") <<endl;

	if(embedding.get("test") != nullptr ){

		cout<<embedding.get("test")<<endl;

	}

	if(embedding.get("dfghtgh") != nullptr ){

		cout<<embedding.get("dfghtgh")<<endl;

	}
	
	



	begin = clock();

	nb_embedded_words_in_voc(embedding , cf);

	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<< "Time to compute nb of words: "<< elapsed_secs <<endl;

	*/
	//launch_experience(lambda , lambda_step , k);

	/*

	vector< vector< pair<int,double> > > results;
	vector< vector<string> > collection;
	vector< vector<string> > queries;
	unordered_map <string,int> cf;
	unordered_map <string,int> df;
	

	clock_t begin = clock();

	//write_all_files("string_content" , "queries.txt" , "cf_file" , "df_file");
	fast_setup("data/collection/string_content" , "data/queries/queries.txt" , "data/stats/cf_file" , "data/stats/df_file" , collection , queries , cf , df);

	clock_t end = clock();
  	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<< "Time to do all the stuff : "<< elapsed_secs <<endl;

	//size_t nb_words = display_stuff(collection , queries , cf , df);

	size_t nb_words = get_size_collection(cf);

	begin = clock();

	//results = basic_language_model(queries , collection , k);
	results = language_model(queries , collection , cf , nb_words , k , lambda);
	
	end = clock();
  	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<< "Time to perform all the queries : "<< elapsed_secs <<endl;
	
	size_t size = 0;

	for(unsigned int i = 0 ; i < results.size() ; i ++){

		size += results[i].size();

	}


	write_res_file(results , "data/res/results" , "CHIC-");
	*/
	
	return 0;

}
