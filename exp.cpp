#include <stdio.h>
#include <stdlib.h> 
#include <omp.h>
#include <iostream>
#include <dirent.h> 
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

//#include <boost/unordered_map.hpp> // unorder_map // check changement // changement 2

using namespace std;

typedef unordered_map<string,int>::const_iterator testit;
typedef unordered_map<int,double>::const_iterator map_iter;


//Check if the input char is a number or a letter 
inline
bool isValidChar(const char c){

	return (( c > 47 && c < 58 ) || ( c > 64 && c < 91 ) || ( c > 96 && c < 123 ));

}

//Check if the input char is a number 
inline
bool isNumber(const char c){

	return (( c > 47 && c < 58 ) );

}

//Check if the input char is a number 
inline
bool isSmallNumber(const char c){

	return (( c > 47 && c < 51 ) );

}

//Check if the input char is a majuscule letter
//inline
//bool isMajChar(const char c){

	//return ( c > 64 && c < 91 );

//}


//Check if the token is valid : no more that 4 digits and no more than 3 identical successiv terms
inline
bool isValidtoken(const string &token){

	if(token=="\0"){return false;}

	if(token.size() > 20){return false;}
	
	int number_count = 0;
	int number_identical_successiv_char = 1;
	
	for(unsigned int i = 0 ; i < token.size() ; i ++){

		if( isNumber(token[i]) ){number_count++;}

		if(i > 0 && token[i]==token[i-1]){number_identical_successiv_char++;}

		if(number_identical_successiv_char > 3 || number_count > 4){return false;}

		if(i > 0 && token[i]!=token[i-1]){number_identical_successiv_char = 1;}

	}

	return true;

}

//Check if the token is valid : no more that 4 digits and no more than 3 identical successiv terms
bool isValidtokendisplay(const string &token){

	int number_count = 0;
	int number_identical_successiv_char = 1;
	
	for(unsigned int i = 0 ; i < token.size() ; i ++){

		if( isNumber(token[i]) ){number_count++;}

		if(i > 0 && token[i]==token[i-1]){number_identical_successiv_char++;}

		if(number_identical_successiv_char > 3 || number_count > 4){return false;}

		if(i > 0 && token[i]!=token[i-1]){number_identical_successiv_char = 1;}

	}

	cout<<"Considered term : "<< token<< " associated number of digits : "<< number_count <<endl;

	return true;

}

//display a string and its associated number of digits 
void nb_digits(const string &token){

	int number_count = 0;

	for(unsigned int i = 0 ; i < token.size() ; i++){

		if( isNumber(token[i]) ){number_count++;}

	}

	cout<<"Considered term : "<< token<< " associated number of digits : "<< number_count <<endl;

}

//Same as before but for each term of a document
void nb_digits(const vector<string> &document){

	for(unsigned int i = 0 ; i < document.size() ; i++){

		nb_digits(document[i]);

	}

}

//Same as before but for the entire collection
void nb_digits(const unordered_map< int , vector<string> > &collection){


	auto iterator = collection.begin();
	
	while(iterator != collection.end()){

		nb_digits(iterator->second);
		iterator++;

	}

}

//Compare the two pairs in input by comparing their second element
bool compare_pairs(const pair<int, double> &p1 , const pair<int, double> &p2) { return (p1.second < p2.second); }


//Displays a document/query
void display(const vector<string> &sentence){

	for(unsigned int i = 0 ; i < sentence.size() ; i++){
		cout<<sentence[i]<<" ";
	}
	cout<<endl;
}


//display the contend of an unordered_map<string,double>
void display_map(const unordered_map<string , double> &map){

	auto iterator = map.begin();

	while(iterator != map.end()){

		cout<< iterator->first << " : " << iterator->second <<endl;
		iterator++;

	}

}



//display the contend of an unordered_map<string,unordered_map<string,double>>
void display_map_map(const unordered_map< string , unordered_map<string,double> > &super_map){

	auto iterator = super_map.begin();

	while(iterator != super_map.end()){

		cout<<"Terme : "<< iterator->first <<endl;
		display_map( iterator->second );

		iterator++;

	}

}


//Aggregate two maps together
void aggregate_map( unordered_map<string,double> &map1 , const  unordered_map<string,double>  &map2){

	auto iterator = map2.begin();

	while(iterator != map2.end()){


		if(map1.find(iterator->first)==map1.end()){map1[iterator->first]=iterator->second;}
		iterator++;

	}

}



//Aggregate two maps together
void aggregate_map( unordered_map< string , unordered_map<string,double> > &map1 , const unordered_map< string , unordered_map<string,double> > &map2){

	auto iterator = map2.begin();

	while(iterator != map2.end()){


		if(map1.find(iterator->first)==map1.end()){map1[iterator->first]=iterator->second;}
		iterator++;

	}

}


//Deletes any big number in the collection
void check_numbers(unordered_map< int , vector<string>  > &collection){
	
	auto iterator = collection.begin();
	
	while(iterator != collection.end()){

		if(iterator->second.size() > 0){		

			for(unsigned int j = 0 ; j < iterator->second.size() ; j++){
				
				if( (iterator->second[j].size() > 4 && isNumber(iterator->second[j][0]) && isNumber(iterator->second[j][4]) )  ){iterator->second.erase(iterator->second.begin() + j);}

			}	

		}
		
		iterator++;

	}

}


//Deletes any empty term of a non empty document in the collection
void check_collection(unordered_map< int , vector<string>  > &collection){

	auto iterator = collection.begin();
	
	while(iterator != collection.end()){

		if(iterator->second.size() > 0){		

			for(unsigned int j = 0 ; j < iterator->second.size() ; j++){

				if(iterator->second[j] == "\0"  ){iterator->second.erase(iterator->second.begin() + j);}			

			}	

		}

		iterator++;

	}

}



//Takes as an input an unordered map of pairs of <int,double> that correspond to the id of the documents and their associated score and returns a vector of pairs containing the k documents that have the highest score
vector< pair<int,double> > kfirst_docs(const unordered_map<int,double> &unsorted, const int k){

	vector< pair<int,double> > res;	

	if(unsorted.size() == 0){return res;}

	pair<int,double> temp;
	
	map_iter p = unsorted.begin();

	while(p != unsorted.end()){

		temp.first = p->first;
		temp.second = p->second;
		res.push_back(temp);
		p++;

	}
	sort(res.begin(), res.end(),compare_pairs);

	if(k != -1 && k < (int)res.size() ){
		res.assign(res.end()-k,res.end());
		
	}

	reverse(res.begin() , res.end());	

	return res;

}


//Returns the terms the have a higher similarity than a given threshold of a given term in the vocabulary 
unordered_map<string,double> closest_terms(const string &term , const unordered_map <string,int> &cf , Embedding &embedding ,  const double &threshold){

	unordered_map<string,double> most_sim;

	if(embedding.get(term.c_str())==nullptr){return most_sim;}

	float cos;
		
	auto iterator = cf.begin();
	while(iterator != cf.end()){
		
		cos =embedding.cosine(term.c_str() , iterator->first.c_str());
		
		if( cos > threshold && term != iterator->first){

			most_sim[iterator->first.c_str()] = cos;

		}
		
		iterator++;

	}
	
	return most_sim;

} 


//Same as before but over the entire vocabulary
unordered_map< string , unordered_map<string,double> > closest_terms(const unordered_map <string,int> &cf , Embedding &embedding ,  const double &threshold){

	unordered_map< string , unordered_map<string,double> > set_most_sim;

	vector<string> temp;

	auto iterator = cf.begin();

	while(iterator != cf.end()){

		temp.push_back(iterator->first);
		iterator++;

	}

	int nthreads, tid;

	#pragma omp parallel private(tid) shared(temp , cf , embedding , threshold)
    {
        tid = omp_get_thread_num();

        unordered_map< string , unordered_map<string,double> > mapLocal;

	

        if (tid == 0)
        {
            nthreads = omp_get_num_threads();
            cout<<"Number of thread: "<<nthreads<<endl;
        }
	
		#pragma omp for schedule(static)
		for(unsigned int i = 0 ; i < temp.size() ; i++){

		//while(iterator != cf.end()){

			//set_most_sim[iterator->first] = closest_terms(iterator->first , cf , embedding , threshold);
			mapLocal[temp[i]] = closest_terms(temp[i] , cf , embedding , threshold);
	
			//iterator++;
		}

		#pragma omp critical
        {
            cout<<"Fusion "<<tid<<" map size="<<mapLocal.size()<<endl;
            set_most_sim.insert(mapLocal.begin(),mapLocal.end());

        }


	}
	
	return set_most_sim;

}



//Same as before but over the queries
unordered_map< string , unordered_map<string,double> > closest_terms(const unordered_map< int , vector<string> > &queries , const unordered_map <string,int> &cf , Embedding &embedding ,  const double &threshold){

	unordered_map< string , unordered_map<string,double> > set_most_sim;

	int compteur = 0;

	auto iterator = queries.begin();

	while(iterator != queries.end()){

		for(unsigned int j = 0 ; j < iterator->second.size() ; j++){

			set_most_sim[iterator->second[j]] = closest_terms(iterator->second[j] , cf , embedding , threshold);

		}

		compteur++;

		cout<<"\rProgress : ("<< compteur <<"/" << queries.size() << ")"<<flush;

		iterator++;

	}

	return set_most_sim;

}


//Returns the terms the have a higher similarity than a given threshold of a given term in the vocabulary 
unordered_map<string,double> closest_terms_sum_query(const vector<string> &query , const unordered_map <string,int> &cf , Embedding &embedding ,  const double &threshold){

	unordered_map<string,double> most_sim;

	if(embedding.check_embedding(query) < 2){return most_sim;}

	float cos;
		
	auto iterator = cf.begin();
	while(iterator != cf.end()){
		
		cos =embedding.cosine(query , iterator->first.c_str());
		
		if( cos > threshold ){

			most_sim[iterator->first.c_str()] = cos;

		}
		
		iterator++;

	}
	
	return most_sim;

} 



//Same as before but over the queries
unordered_map< int , unordered_map<string,double> > closest_terms_sum_query(const unordered_map< int , vector<string> > &queries , const unordered_map <string,int> &cf , Embedding &embedding ,  const double &threshold){

	unordered_map< int , unordered_map<string,double> > set_most_sim;

	int compteur = 0;

	auto iterator = queries.begin();

	while(iterator != queries.end()){


		set_most_sim[iterator->first] = closest_terms_sum_query(iterator->second , cf , embedding , threshold);

		compteur++;

		cout<<"\rProgress : ("<< compteur <<"/" << queries.size() << ")"<<flush;

		iterator++;

	}

	return set_most_sim;

}




//Write a <string,double> map in a file
void write_map(const unordered_map< string , double > &map , const string &file_name){

	ofstream myfile;
  	myfile.open (file_name.c_str());
	
	auto iterator = map.begin();

	while( iterator != map.end() ){

		myfile << iterator->first + " " + to_string(iterator->second) + "\n";
	
		iterator++;

	}

  	myfile.close();

}


//Take a string as input and return a pair of string and double using the delimiter
pair<string,double> read_pair(const string &str, const char delimiter) {

  pair<string,double> res;
  stringstream ss(str); // Turn the string into a stream.
  string tok;
  tok[0] = 0;
  getline(ss, tok, delimiter);
  res.first = tok;
  getline(ss, tok, delimiter);
  res.second = atof(tok.c_str());		
  
  return res;
}


//Read a file containing the map
unordered_map< string , double > read_cos_map_file(const string &file_name , const double &threshold){

	FILE* f = fopen(file_name.c_str(),"r");

	size_t len = 10000;

	char* line = (char*)malloc(len);

	size_t nread;

	unordered_map< string , double > res;

	pair<string,double> temp;

	while (int(nread = getline(&line, &len, f)) != -1) {

		temp = read_pair( string(line) , ' ');
		if(temp.second > threshold){res[temp.first] = temp.second;}
		
	}

	free(line);	

	return res;

}




//Write in a file the set of the closest words of the vocabulary
void write_map_map(const unordered_map< string , unordered_map<string,double> > &set_closest_words , const string &file_name){

	ofstream myfile;
  	myfile.open (file_name.c_str());
	
	auto iterator = set_closest_words.begin();
	unordered_map<string,double>::const_iterator iterator2;

	while( iterator != set_closest_words.end() ){

		myfile << iterator->first + "\n";

		iterator2 = iterator->second.begin();

		while(iterator2 != iterator->second.end()){

			myfile << iterator2->first + " " + to_string(iterator2->second) + " ";
			iterator2++;

		}
		
		myfile << "\n";
	
		iterator++;

	}

  	myfile.close();

}






//Write in a file the set of the closest words of the vocabulary
void write_map_map(const unordered_map< int , unordered_map<string,double> > &set_closest_words , const string &file_name){

	ofstream myfile;
  	myfile.open (file_name.c_str());
	
	auto iterator = set_closest_words.begin();
	unordered_map<string,double>::const_iterator iterator2;

	while( iterator != set_closest_words.end() ){

		myfile << to_string(iterator->first) + "\n";

		iterator2 = iterator->second.begin();

		while(iterator2 != iterator->second.end()){

			myfile << iterator2->first + " " + to_string(iterator2->second) + " ";
			iterator2++;

		}
		
		myfile << "\n";
	
		iterator++;

	}

  	myfile.close();

}




//Takes a string that corresponds to one line of the mapmap file as an input and returns an unordered_map
unordered_map<string,double> read_map(const string &line , const double &threshold){

	unordered_map<string,double> res;
	
	stringstream ss(line);
  	string tok;
	string tok2;
	while(getline(ss, tok, ' ')){

		getline(ss, tok2, ' ');

		if(tok!="\n" && atof(tok2.c_str()) > threshold ){
			
			res[tok] = atof(tok2.c_str());
		}

	}

	return res;

}

//Read a file containing the mapmap
unordered_map< string , unordered_map<string,double> > read_cos_map_map_file(const string &file_name , const double &threshold){

	FILE* f = fopen(file_name.c_str(),"r");

	size_t len = 10000;

	char* line = (char*)malloc(len);
	char* line2 = (char*)malloc(len);

	size_t nread;

	unordered_map< string , unordered_map<string,double> > res;

	while (int(nread = getline(&line, &len, f)) != -1) {

		getline(&line2, &len, f);
		
		res[string(line).erase(string(line).size()-1)] = read_map(line2 , threshold);

	}

	free(line);	
	free(line2);
	
	return res;

}


//Read a file containing the mapmap
unordered_map< int , unordered_map<string,double> > read_cos_sum_query_map_map_file(const string &file_name , const double &threshold){

	FILE* f = fopen(file_name.c_str(),"r");

	size_t len = 10000;

	char* line = (char*)malloc(len);
	char* line2 = (char*)malloc(len);

	size_t nread;

	unordered_map< int , unordered_map<string,double> > res;

	while (int(nread = getline(&line, &len, f)) != -1) {

		getline(&line2, &len, f);
		
		res[atoi( string(line).erase(string(line).size()-1).c_str() )] = read_map(line2 , threshold);

	}

	free(line);	
	free(line2);
	
	return res;

}




//Return the cosine similarity between two terms ; return 0 if one of the term is not in the input map 
double fast_cos(const string &term1 , const string &term2 , const unordered_map< string , unordered_map<string,double> > &set_cos){

	if(term1==term2){return 1;}

	auto iterator = set_cos.find( term1 );

	if(iterator != set_cos.end()){

		auto iterator2 = iterator->second.find(term2);

		if(iterator2 == iterator->second.end()){return 0;}
		
		else{return iterator2->second;}

	}

	else return 0;

}



//Return the neighbouring terms of the input term 
vector<string> neighbors(const vector<string> &query , const unsigned int term_position ,const int size_neighborhood){

	unsigned int pos = 0; 

	if(term_position - size_neighborhood > 0){pos = term_position - size_neighborhood;}

	vector<string> neighborhood;

	while(pos < term_position + size_neighborhood + 1 && pos < query.size()){

		neighborhood.push_back(query[pos]);
		pos++;

	}

	return neighborhood;

}



//Reads a qrel file, deletes the lines corresponding to the non relevant file ans save the new qrel file 
int correct_qrel_file(const string &filenamein , const string &filenameout){
	
	string line;

    ifstream in(filenamein);
    if( !in.is_open())
    {
          cout << "Input file failed to open\n";
          return 1;
    }

    ofstream out(filenameout);

    while( getline(in,line) )
    {
        if(line[line.size()-1] != '0')
            out << line << "\n";
    }
    in.close();
    out.close();    
    return 0;

}







///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//Write the collection/set of queries in input in a file
void write_xml_collection_file(const unordered_map< int , vector<string> > &collection, const string &file_name){

	
	ofstream myfile;
  	myfile.open (file_name.c_str());

	myfile << "<matrix>\n";
	
	auto iterator = collection.begin();

	while( iterator != collection.end() ){

		myfile << "<DOC>\n";
		myfile << "<DOCNO>"<<iterator->first<<"</DOCNO>\n";
		myfile << "<TEXT>";
		if(iterator->second.size() != 0){
		
			for(unsigned int j = 0 ; j < iterator->second.size()-1 ; j++){
			
				myfile << iterator->second[j] << " ";

			}

			myfile << iterator->second[ iterator->second.size()-1 ];

		}
		
		myfile << "</TEXT>\n";
		myfile << "</DOC>\n";					

		iterator++;

	}	

	myfile << "</matrix>\n";
	myfile.close();

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




//Write the unordered map of <string,int> in input in a file 
void write_tf_file(const unordered_map <string,int>  &cf , const string &file_name){

	ofstream myfile;
  	myfile.open (file_name.c_str());

	testit p = cf.begin();

	while( p!=  cf.end() ){

		myfile << p->first << " " << p->second << "\n";
		p++;

	}
	
    myfile.close();

}


//Take a string as input and return a vector of strings using the delimiter
vector<string> easy_split(const string &str, const char delimiter) {

  vector<string> internal;
  stringstream ss(str); // Turn the string into a stream.
  string tok;
  tok[0] = 0;
  while(getline(ss, tok, delimiter)) {
			internal.push_back(tok);		
  }
  
  return internal;
}

//Take a string as input and return a pair<string,int>
pair<string,int> tf_split(const string &str, const char delimiter) {

  pair<string,int> internal;
  stringstream ss(str); // Turn the string into a stream.
  string tok;
  tok[0] = 0;
  getline(ss, tok, delimiter);
  internal.first = tok;		
  getline(ss, tok, delimiter);
  internal.second = stoi(tok);
  
  return internal;
}





//Take a string as input and return a vector of strings using the delimiter also get rides of the caracters that are not numbers and letters 
vector<string> split(const string &str, const char delimiter) {

  vector<string> internal;
  stringstream ss(str); // Turn the string into a stream.
  string tok;
  tok[0] = 0;
  unsigned int i = 0;
  while(getline(ss, tok, delimiter)) {
		while(tok.size() > 1 && !isValidChar(tok[0]) ){
			tok.erase(0);
		}
		if(tok.size()>1){
			i = 0;
			while(i < tok.size()){
				//if(!isValidChar(tok[i]) && tok[i]!='-'){
				if( !isValidChar(tok[i]) ){
					internal.push_back(tok.substr(0,i));
					tok = tok.substr(i+1);
					while(tok.size() > 1 && !isValidChar(tok[0]) ){
						tok.erase(0);
					}
					i = -1;
				}
				i++;
			}
		}

		while(tok.size() > 1 && !isValidChar(tok[0]) ){
			tok.erase(0);
		}

		while(tok.size() > 1 && !isValidChar(tok[tok.size()-1])){
			tok.erase(tok.size()-1);
		}
		if(tok.size() != 0 && isValidChar(tok[0])){
			internal.push_back(tok);		
		}
  }
  
  return internal;
}








//Take a string as input and return a vector of strings using the delimiter also get rides of the caracters that are not numbers and letters and changes majuscules into minuscules
vector<string> split_maj(const string &str, const char delimiter) {

  vector<string> internal;
  stringstream ss(str); // Turn the string into a stream.
  string tok;
  string sub_tok;
  tok[0] = 0;
  unsigned int i = 0;
  while(getline(ss, tok, delimiter)) {
		while(tok.size() > 1 && !isValidChar(tok[0]) ){
			tok.erase(0);
		}
		if(isMajChar(tok[0])){tok[0] = tok[0] + 32;}
		if(tok.size()>1){
			i = 0;
			while(i < tok.size()){
				if(isMajChar(tok[i])){tok[i] = tok[i] + 32;}
				//if(!isValidChar(tok[i]) && tok[i]!='-'){
				if( !isValidChar(tok[i]) ){
					sub_tok = tok.substr(0,i);
					if(isValidtoken(sub_tok)){internal.push_back(sub_tok);}
					tok = tok.substr(i+1);
					while(tok.size() > 1 && !isValidChar(tok[0]) ){
						tok.erase(0);
					}
					i = -1;
				}
				i++;
			}
		}

		while(tok.size() > 1 && !isValidChar(tok[0]) ){
			tok.erase(0);
		}

		while(tok.size() > 1 && !isValidChar(tok[tok.size()-1])){
			tok.erase(tok.size()-1);
		}
		if(tok.size() != 0 && isValidChar(tok[0]) && isValidtoken(tok) ){
			internal.push_back(tok);		
		}
  }
  
  return internal;
}



//Take a string as input and return a vector of strings using the delimiter also get rides of the caracters that are not numbers and letters and changes majuscules into minuscules
vector<string> skip(const string &str, const char delimiter) {

  vector<string> internal;
  stringstream ss(str); // Turn the string into a stream.
  string tok;
  string sub_tok;
  tok[0] = 0;
  unsigned int i = 0;
  while(getline(ss, tok, delimiter)) {
		i++;
  }
  
  return internal;
}


//Take a string as input and return a vector of strings using the delimiter also get rides of the caracters that are not numbers and letters and changes majuscules into minuscules
vector<string> split_maj_display(const string &str, const char delimiter) {

  vector<string> internal;
  stringstream ss(str); // Turn the string into a stream.
  string tok;
  tok[0] = 0;
  unsigned int i = 0;
  while(getline(ss, tok, delimiter)) {
		while(tok.size() > 1 && !isValidChar(tok[0]) ){
			tok.erase(0);
		}
		if(isMajChar(tok[0])){tok[0] = tok[0] + 32;}
		if(tok.size()>1){
			i = 0;
			while(i < tok.size()){
				if(isMajChar(tok[i])){tok[i] = tok[i] + 32;}
				//if(!isValidChar(tok[i]) && tok[i]!='-'){
				if( !isValidChar(tok[i]) &&isValidtokendisplay(tok.substr(0,i)) ){

					internal.push_back(tok.substr(0,i));
					tok = tok.substr(i+1);
					while(tok.size() > 1 && !isValidChar(tok[0]) ){
						tok.erase(0);
					}
					i = -1;
				}
				i++;
			}
		}

		while(tok.size() > 1 && !isValidChar(tok[0]) ){
			tok.erase(0);
		}

		while(tok.size() > 1 && !isValidChar(tok[tok.size()-1])){
			tok.erase(tok.size()-1);
		}
		if(tok.size() != 0 && isValidChar(tok[0]) && isValidtokendisplay(tok) ){
			internal.push_back(tok);		
		}
  }
  
  return internal;
}


//Reads a collection/set of queries in the input file
unordered_map< int , vector<string> > read_some_file(const string &file_name){

  FILE* f = fopen(file_name.c_str(),"r");

	size_t len = 10000;

	char* line = (char*)malloc(len);

	size_t nread;

	int current_iteration = 0;

	unordered_map< int , vector<string> > collection;


	while (int(nread = getline(&line, &len, f)) != -1) {

		if( (current_iteration)==377000){collection[current_iteration] =  split_maj(string(line) , ' ');}
		else{ skip(string(line) , ' ') ;}
		current_iteration++;
		
	}

	free(line);	
	
	//check_numbers(collection);
	check_collection(collection);

	auto iterator = collection.begin();

	//for(unsigned int i = 0 ; i < collection.size() ; i++){
		
	while( iterator != collection.end() ){	
		cout<<"Document "<<iterator->first<<" :"<<endl;
		display(iterator->second);
		cout<<endl<<endl;

	}

	return collection;
	
}


//Reads a collection/set of queries in the input file
unordered_map< int , vector<string> > read_file(const string &file_name){

  FILE* f = fopen(file_name.c_str(),"r");

	size_t len = 10000;

	char* line = (char*)malloc(len);

	size_t nread;

	unordered_map< int , vector<string> > collection;

	int current_iteration = 0;

	while (int(nread = getline(&line, &len, f)) != -1) {

		collection[current_iteration] =  split_maj(string(line) , ' ');

		current_iteration++;
		
	}

	free(line);	

	
	//check_numbers(collection);
	check_collection(collection);

	return collection;
	
}


//Same as before but using the function easy_split on the lines of the file
unordered_map< int , vector<string> > read_clean_file(const string &file_name){

  FILE* f = fopen(file_name.c_str(),"r");

	size_t len = 10000;

	char* line = (char*)malloc(len);

	size_t nread;

	unordered_map< int , vector<string> > collection;

	int current_iteration = 0;

	while (int(nread = getline(&line, &len, f)) != -1) {

		collection[current_iteration] = easy_split(string(line) , ' ');

		current_iteration++;

	}

	free(line);

	return collection;
	
}


//Same as before but using the function tf_split on the lines of the file
unordered_map <string,int> read_tf_file(const string &file_name){


	 FILE* f = fopen(file_name.c_str(),"r");

	size_t len = 10000;

	char* line = (char*)malloc(len);

	size_t nread;

	unordered_map <string,int> cf;
	

	while (int(nread = getline(&line, &len, f)) != -1) {


		cf.insert( tf_split(string(line) , ' ') );


	}
	
	free(line);

    return cf;

}



//Builds the vocabulary given the collection in input and stores it in an unordered_map of <string,int> the string being the words and the int their associated collection frequency
unordered_map <string,int> build_cf(const unordered_map< int , vector<string> > &collection){

	unordered_map <string,int>  cf;

	auto iterator = collection.begin();

	while( iterator != collection.end() ){

		for(unsigned int j = 0 ; j < iterator->second.size() ; j++){
					
			if(iterator->second[j].size() > 1 || isValidChar(iterator->second[j][0])){
				cf[ iterator->second[j] ]++;
			}
		}

		iterator++;
	
	}

	return cf;
	
}

//Given the collection and the vocabulary in input computes and store the document frequency of the terms of the vocabulary in an unordered_map
unordered_map <string,int> build_df(const unordered_map< int , vector<string> > &collection, const unordered_map <string,int> &cf){

	unordered_map <string,int>  df;

	testit p = cf.begin();
	
	pair<string,int> temp("",0);

	while( p!=  cf.end() ){
		
		temp.first = p->first;
		df.insert(temp);
		p++;

	}

	vector<string> doc_memory;

	auto iterator = collection.begin();

	while( iterator != collection.end() ){

		for(unsigned int j = 0 ; j < iterator->second.size() ; j++){

			if( count( doc_memory.begin() , doc_memory.end() , iterator->second[j] ) == 0 && (df.find(iterator->second[j]) != df.end())){

				df[ iterator->second[j] ]++;
				doc_memory.push_back(iterator->second[j]);

			}

		}
		
		doc_memory.clear();

		iterator++;

	}

	return df;
	
}





//Reads a qrel file, separates the results between a training and a test set
int separate_qrel_file(const string &filename , const string &filename_train , const string &filename_test , const unordered_map<int , vector<string> > &collection){

	string line;

    ifstream in(filename);
    if( !in.is_open() )
    {
          cout << "Input file failed to open\n";
          return 1;
    }

    ofstream train(filename_train);

	ofstream test(filename_test);

    while( getline(in,line) )
    {

        if( collection.find(stoi( easy_split(line, ' ')[2] )) != collection.end()  ){

			train << line << "\n";	
		
		}

		else{test << line << "\n";}
            
    }
    in.close();
    train.close();
	test.close();    
    return 0;

}


//Return the tf of a term in a document
int term_freq(const vector<string> &doc , const string &term){

	return count(doc.begin() , doc.end() , term);

}

//Return the collection frequency corresponding to the term in input
int coll_freq(const unordered_map <string,int>  &cf , const string &term){

	testit q = cf.find(term); 

	if(q == cf.end() || q->second == -1){return 0;}
	else{return q->second;}

}

//Return the probability that the query was generated by the model of the document using only term frequency in the document (no smoothing)
double basic_language_model(const vector<string> &query , const vector<string> &document){
	
	size_t doc_length = document.size();

	if(doc_length == 0){return 0;}

	double proba = 0;

	vector<string> check;

	for(unsigned int i = 0 ; i < query.size() ; i++){
	
		if( find(check.begin() , check.end() , query[i]) == check.end() ){
			proba += log( 1 + (double)term_freq(document , query[i])/doc_length );
			check.push_back(query[i]);

		}

	}

	return proba;

}

//Return a list of the probability that the query was generated by the model of each of the documents using only term frequency in the document (no smoothing)
unordered_map <int,double> basic_language_model(const vector<string> &query , const unordered_map < int , vector<string> > &collection){
	
	unordered_map <int,double> list_doc;

	double proba;

	auto iterator = collection.begin();

	while(iterator != collection.end() ){

		proba = basic_language_model(query , iterator->second);
		if(proba > 0){list_doc[iterator->first]=proba;}

		iterator++;

	}

	return list_doc;

}


//Same as before but with all the queries and sort documents by their score
vector< vector< pair<int,double> > > basic_language_model(const unordered_map <int , vector<string> > &queries , const unordered_map < int , vector<string> > &collection , const int k){

	vector< vector< pair<int,double> > > list_docs(queries.size());
	
	auto iterator = queries.begin();	

	while( iterator != queries.end() ){

		//cout<<"\rProcessing query "<<i+1<<"/"<<queries.size()<<flush;

		list_docs[iterator->first] = kfirst_docs( basic_language_model(iterator->second , collection) , k);	

		iterator++;

	}

	return list_docs;

}


//Same as before but with a smoothing that takes into account the collection frequency
double Hiemstra_language_model(const vector<string> &query , const vector<string> &document , const unordered_map <string,int>  &cf , const int collection_size , const double lambda){
	
	size_t doc_length = document.size();

	if(doc_length == 0){return 0;}

	double doc_proba;

	double coll_proba;

	double proba = 0;

	for(unsigned int i  = 0 ; i < document.size() ; i++){

		if(document[i]=="\0"){cout<<"|"<<document[i]<<"|"<<" and size of the doc : "<<document.size()<<endl;}

	}

	for(unsigned int i = 0 ; i < query.size() ; i++){

			coll_proba =  (1 - lambda)*( (double)coll_freq(cf , query[i])/collection_size );

			if(coll_proba!=0){
				
				doc_proba = (lambda)*( (double)term_freq(document , query[i])/doc_length );
				
				proba += log(1 + doc_proba/coll_proba)/log(2);

			}
	
	}

	return proba;

}



//Same as before but with a smoothing that takes into account the collection frequency
unordered_map <int,double> Hiemstra_language_model(const vector<string> &query , const unordered_map < int , vector<string> > &collection , const unordered_map <string,int>  &cf , const int collection_size , const double lambda){
	
	unordered_map <int,double> list_doc;

	double proba;

	auto iterator = collection.begin();

	while(iterator != collection.end()){

		proba = Hiemstra_language_model(query , iterator->second , cf , collection_size , lambda);
		if(proba!=0){list_doc[iterator->first]=proba;}

		iterator++;
		
	}

	return list_doc;

}



//Same as before but with all the queries and sort documents by their score
vector< vector< pair<int,double> > > Hiemstra_language_model(const unordered_map< int , vector<string> > &queries , const unordered_map < int , vector<string> > &collection , const unordered_map <string,int>  &cf , const int collection_size , const int k , const double lambda){

	vector< vector< pair<int,double> > > list_docs(queries.size());

	auto iterator = queries.begin();
	
	while(iterator != queries.end()){

		//cout<<"\rProcessing query "<<i+1<<"/"<<queries.size()<<flush;

		list_docs[iterator->first] = kfirst_docs( Hiemstra_language_model(iterator->second , collection , cf , collection_size , lambda) , k );	

		iterator++;

	}

	cout<<endl;

	return list_docs;

}




//Computes the log of the probability that query was generated by the document w.r.t dirichlet language model pdir(query | document)
inline
double Dirichlet_language_model(const double &mu , const vector<string> &query , const vector<string> &document , const unordered_map <string,int>  &cf , const int collection_size){

	if(document.size()==0){return 0;}

	double res = 0;

	//double collection_proba;

	for(unsigned int i = 0 ; i < query.size() ; i++){

		//collection_proba = (double)coll_freq(cf , query[i])/collection_size;	

		//if(coll_freq(cf , query[i]) != 0){

			//res += log(mu/(mu+document.size()));

			//if(term_freq(document , query[i]) != 0){	
		
				//res += log( 1 + term_freq(document , query[i])/(mu*(double)coll_freq(cf , query[i])/collection_size)) + log(mu/(mu+document.size()));
		
			//}
					
		//}



		if(coll_freq(cf , query[i]) != 0){

			res += log( ( term_freq(document , query[i]) + mu*((double)coll_freq(cf , query[i])/collection_size))/(document.size() + mu) );
		
		}

	}


	return res;

}



//Same as before but over the entire collection
unordered_map <int,double> Dirichlet_language_model(const double &mu , const vector<string> &query , const unordered_map< int , vector<string> > &collection , const unordered_map <string,int>  &cf , const int collection_size ){
	
	unordered_map <int,double> list_doc;

	double proba;

	auto iterator = collection.begin();

	while(iterator != collection.end()){

		proba = Dirichlet_language_model(mu , query , iterator->second , cf , collection_size);
		if(proba!=0){list_doc[iterator->first]=proba;}

		iterator++;
		
	}

	return list_doc;

}
 

//Same as before but with all the queries and sort documents by their score
vector< vector< pair<int,double> > > Dirichlet_language_model(const double &mu , const unordered_map< int , vector<string> > &queries , const unordered_map < int , vector<string> > &collection , const unordered_map <string,int>  &cf , const int k , const int collection_size){

	vector< vector< pair<int,double> > > list_docs(queries.size());
	
	auto iterator = queries.begin();
	
	while(iterator != queries.end()){

		//cout<<"\rProcessing query "<<i+1<<"/"<<queries.size()<<flush;

		list_docs[iterator->first] =  kfirst_docs( Dirichlet_language_model(mu , iterator->second , collection , cf , collection_size) , k );

		iterator++;

	}

	return list_docs;

}






//Read the sum of the similarities between a term and all the terms of a given set of terms by using a map
double fast_sum_proba_pairs(const string &term , const vector<string> &set , const unordered_map< string , unordered_map<string,double> > &cosine_map){

	double sum = 0;

	for(unsigned int i = 0 ; i < set.size() ; i++){

		sum += fast_cos(term , set[i] , cosine_map);

	}

	return sum;

}



//Read a map and computes the sum of all the cosine similarities of a given term 
double fast_sum_cos_sim(const string &term , const unordered_map< string , unordered_map<string,double> > &cosine_map){

	auto iterator = cosine_map.find( term );

	if(iterator != cosine_map.end()){

		auto iterator2 = iterator->second.begin();

		double sum = 0;

		while( iterator2 != iterator->second.end() ){

			sum+=iterator2->second;
			iterator2++;

		}

		return sum;

	}

	else{return 0;}

}


//Read a map and computes and store into a map the sum of all the cosine similarities of all the terms in the input map
unordered_map<string , double> all_fast_sum_cos_sim(const unordered_map< string , unordered_map<string,double> > &cosine_map){

	unordered_map<string , double> res;
	
	auto iterator = cosine_map.begin();

	while(iterator != cosine_map.end()){

		
		res[iterator->first] = fast_sum_cos_sim( iterator->first , cosine_map );
		iterator++;

	} 

	return res;

} 



double fast_sum_translation_proba(const string &term1 , const string &term2 , const unordered_map< string , unordered_map<string,double> > &sum_proba ){

	auto iterator = sum_proba.find(term1);

	if(iterator == sum_proba.end()){return 0;}

	auto iterator2 = iterator->second.find(term1);

	if(iterator2 == iterator->second.end() ){

		return 0;

	}

	else{return iterator2->second;}

}


//Computes the translation probability of a term2 into term1 : p( term1 | term2 )
double translation_proba(const string &term1 , const string &term2 , const unordered_map<string , double> &sum_cosine_map , const unordered_map< string , unordered_map<string,double> > &cosine_map ){

	//if(term1==term2){return 1;}

	//else{return 0;}

	auto iterator = sum_cosine_map.find(term2);

	if(iterator == sum_cosine_map.end() || iterator->second == 0){
		
		return fast_cos(term1 , term2 , cosine_map);

	}

	else{

		return  fast_cos(term1 , term2 , cosine_map)/iterator->second;	

	}

}




//Computes the probability that term was generated by the document p(term | document)
double proba_doc_generate_term(const string &term , const vector<string> &document , const unordered_map<string , double> &sum_cosine_map , const unordered_map< string , unordered_map<string,double> > &cosine_map ){
	
	double proba = 0;

	for(unsigned int i = 0 ; i < document.size() ; i++){

		proba+=translation_proba(term , document[i] , sum_cosine_map , cosine_map )/document.size();

	}

	return proba;

}



//Computes the probability that term was generated by the document p(term | document) by taking into account the self translation probability
double proba_doc_generate_term(const string &term , const vector<string> &document , const unordered_map<string , double> &sum_cosine_map , const unordered_map< string , unordered_map<string,double> > &cosine_map , const double &alpha ){
	
	double proba = 0;

	//double self;

	for(unsigned int i = 0 ; i < document.size() ; i++){

		if(term == document[i]){

			proba+= (alpha + (1 - alpha)*translation_proba(term , document[i] , sum_cosine_map , cosine_map ))/document.size();
			
			//proba+= alpha/document.size();

		}

		else{

			proba+=(1 - alpha)*translation_proba(term , document[i] , sum_cosine_map , cosine_map )/document.size();

			//self = translation_proba(document[i] , document[i] , sum_cosine_map , cosine_map );

			//if(self==1){proba += (1 - alpha)*translation_proba(term , document[i] , sum_cosine_map , cosine_map )/(document.size());}
	
			//else{proba += (1 - alpha)*translation_proba(term , document[i] , sum_cosine_map , cosine_map )/((1-self)*document.size());}

		}
		
	}

	return proba;

}



//Computes the log of the probability that query was generated by the document w.r.t dirichlet language model  pdir(query | document)
double Dirichlet_embedding_model(const double &mu , const vector<string> &query , const vector<string> &document , const unordered_map <string,int>  &cf , const unordered_map<string , double> &sum_cosine_map , const unordered_map< string , unordered_map<string,double> > &cosine_map , const int collection_size , const double &alpha ){

	double res = 0;

	double proba;

	int collection_freq;

	for(unsigned int i = 0 ; i < query.size() ; i++){

		//if(term_freq(document , query[i]) != 0){	

			//res += log( 1 + ( proba_doc_generate_term(query[i] , document , sum_cosine_map , cosine_map)*document.size() )/(mu*(double)coll_freq(cf , query[i])/collection_size) ) + log(mu/(mu+document.size()));
		
		//}
	
		proba = proba_doc_generate_term(query[i] , document , sum_cosine_map , cosine_map , alpha);

		collection_freq = coll_freq(cf , query[i]);

		if(collection_freq != 0 && proba != 0){

			res += log( ( (proba*document.size()) + mu*((double)collection_freq/collection_size))/(document.size() + mu) );
		
		}

		else if(collection_freq != 0 && proba == 0){

			res += log( (double)collection_freq/collection_size );

		}

		else if(collection_freq == 0 && proba != 0){

			res += log( proba );

		}

	}

	return res; 

}



//Same as before but over the entire collection
unordered_map <int,double> Dirichlet_embedding_model(const double &mu , const vector<string> &query , const unordered_map< int , vector<string> > &collection , const unordered_map <string,int>  &cf , const unordered_map<string , double> &sum_cosine_map , const unordered_map< string , unordered_map<string,double> > &cosine_map , const int collection_size , const double &alpha){
	
	unordered_map <int,double> list_doc;

	double proba;

	auto iterator = collection.begin();

	while( iterator != collection.end() ){

		proba = Dirichlet_embedding_model(mu , query , iterator->second , cf , sum_cosine_map , cosine_map , collection_size , alpha );
		if(proba!=0){list_doc[iterator->first]=proba;}
		
		iterator++;

	}

	return list_doc;

}


//Same as before but with all the queries and sort documents by their score
vector< vector< pair<int,double> > > Dirichlet_embedding_model(const double &mu , const unordered_map< int , vector<string> > &queries , const unordered_map< int , vector<string> > &collection , const unordered_map <string,int>  &cf , const unordered_map<string , double> &sum_cosine_map , const unordered_map< string , unordered_map<string,double> > &cosine_map , const int k , const int collection_size , const double &alpha ){

	vector< vector< pair<int,double> > > list_docs(queries.size());

	auto iterator = queries.begin();
	
	while(iterator != queries.end() ){

		list_docs[iterator->first] = kfirst_docs( Dirichlet_embedding_model(mu , iterator->second , collection , cf , sum_cosine_map , cosine_map , collection_size , alpha ) , k );

		iterator++;

	}

	cout<<endl;

	return list_docs;

}


/*


//Computes the translation probability of a term2 into term1 : p( term1 | term2 )
double sum_translation_proba(const string &term1 , const string &term2 , const unordered_map<string , double> &sum_cosine_map , const unordered_map< string , unordered_map<string,double> > &cosine_map ){


	if(term1==term2){return 0;}

	auto iterator = cosine_map.find(term2);

	int input;

	if(iterator == cosine_map.end()){return 0;}
	
	double sum = -1.0*translation_proba(term1 , term2 , sum_cosine_map , cosine_map );

	if(sum == 0){return 0;}

	auto iterator2 = iterator->second.begin();

	if(iterator->second.size() == 1){return 1;}

	cout<<"Size of the sum : "<< iterator->second.size() <<endl;

	cout<<"First value of the sum : "<< sum <<endl;

	while(iterator2 != iterator->second.end()){

		sum += translation_proba(iterator2->first , term2 , sum_cosine_map , cosine_map );

		iterator2++;
	
	}

	cout<<"Sum translation proba between "<< term1 <<" and "<< term2 << " : "<< sum <<endl;

	cin>>input;

	return sum;

}


//Computes the probability that term was generated by the document p(term | document) by taking into account the self translation probability
unordered_map<string , double> save_sum_proba(const string &term , const unordered_map<string,int> &cf , const unordered_map<string , double> &sum_cosine_map , const unordered_map< string , unordered_map<string,double> > &cosine_map ){
	
	unordered_map<string , double > res;

	double sum;

	auto iterator = cf.begin();
	
	while(iterator != cf.end()){
		
		sum = sum_translation_proba(term , iterator->first , sum_cosine_map , cosine_map );

		if(sum != 0){

			res[iterator->first] = sum;

		}

		iterator++;
	
	}

	return res;

}


//Same as before but with all the queries and sort documents by their score
unordered_map< string , unordered_map<string , double > > save_sum_proba(const unordered_map< int , vector<string> > &queries , const unordered_map <string,int> &cf ,  const unordered_map<string , double> &sum_cosine_map , const unordered_map< string , unordered_map<string,double> > &cosine_map ){

	unordered_map< string , unordered_map<string , double > > res;

	vector<string> temp;

	//auto iterator = queries.begin();

	int nthreads, tid;

	#pragma omp parallel private(tid)  shared(queries , cf , sum_cosine_map , cosine_map) num_threads(1)
	{
		tid = omp_get_thread_num();
		unordered_map< string , unordered_map<string , double > > mapLocal;

		if (tid == 0)
		{
		    nthreads = omp_get_num_threads();
		    cout<<"Number of thread: "<<nthreads<<endl;
		}

		#pragma omp for schedule(static)
		for(unsigned int j = 0 ; j < queries.size() ; j++){

		//while(iterator != queries.end()){

			//cout<<"\rProcessing query "<<iterator->first<<"/"<<queries.size()<<flush;
	
			//for(unsigned int i = 0 ; i < iterator->second.size() ; i++){

				//res[iterator->second[i]] = save_sum_proba(iterator->second[i] , cf , sum_cosine_map , cosine_map );
					
			for(unsigned int i = 0 ; i < queries.find(j)->second.size() ; i++){
		
				//res[queries.find(j)->second[i]] = save_sum_proba(queries.find(j)->second[i] , cf , sum_cosine_map , cosine_map );

				mapLocal[queries.find(j)->second[i]] = save_sum_proba(queries.find(j)->second[i] , cf , sum_cosine_map , cosine_map );

			}

		}

		#pragma omp critical
        {
            cout<<"Fusion "<<tid<<" map size="<<mapLocal.size()<<endl;
            res.insert(mapLocal.begin(),mapLocal.end());

        }

	}
	//iterator++;

	return res;	

}

*/

//Read a map and computes the sum of all the cosine similarities of a given term 
double fast_cos_sum_queries(const int query_id , const string &term , const unordered_map< int , unordered_map<string,double> > &all_cos_sum_queries){

	auto iterator = all_cos_sum_queries.find( query_id );

	if(iterator != all_cos_sum_queries.end()){

		auto iterator2 = iterator->second.find(term);

		if(iterator2 != iterator->second.end() ){

			return iterator2->second;

		}

	}

	return 0;

}



//Computes the translation probability of a term2 into phrase : p( phrase | term2 )
double translation_proba(const int query_id , const string &term2 , const unordered_map<string , double> &sum_cosine_map , unordered_map< int , unordered_map<string,double> > &all_cos_sum_queries){

	//if(term1==term2){return 1;}

	//else{return 0;}

	auto iterator = sum_cosine_map.find(term2);

	if(iterator == sum_cosine_map.end() || iterator->second == 0){
		
		return fast_cos_sum_queries(query_id , term2 , all_cos_sum_queries);

	}

	else{

		return  fast_cos_sum_queries(query_id , term2 , all_cos_sum_queries)/iterator->second;	

	}

}




//Computes the probability that the phrase was generated by the document p(term | document)
double proba_doc_generate_term(const int query_id , const vector<string> &document , const unordered_map<string , double> &sum_cosine_map , unordered_map< int , unordered_map<string,double> > &all_cos_sum_queries){
	
	double proba = 0;

	for(unsigned int i = 0 ; i < document.size() ; i++){

		proba+=translation_proba( query_id , document[i] , sum_cosine_map , all_cos_sum_queries )/document.size();

	}

	return proba;

}



//Computes the log of the probability that query was generated by the document w.r.t dirichlet language model  pdir(query | document)
double Dirichlet_embedding_model_Qexp(const double &mu , const vector<string> &query , const int query_id , const vector<string> &document , const unordered_map <string,int>  &cf , const unordered_map<string , double> &sum_cosine_map , const unordered_map< string , unordered_map<string,double> > &cosine_map , unordered_map< int , unordered_map<string,double> > &all_cos_sum_queries , const int collection_size){

	double res = 0;

	double proba;

	int collection_freq;

	for(unsigned int i = 0 ; i < query.size() ; i++){
	
		proba = proba_doc_generate_term(query[i] , document , sum_cosine_map , cosine_map);

		collection_freq = coll_freq(cf , query[i]);

		if(collection_freq != 0 && proba != 0){

			res += log( ( (proba*document.size()) + mu*((double)collection_freq/collection_size))/(document.size() + mu) );
		
		}

		else if(collection_freq != 0 && proba == 0){

			res += log( (double)collection_freq/collection_size );

		}

		else if(collection_freq == 0 && proba != 0){

			res += log( proba );

		}

	}

	proba = proba_doc_generate_term( query_id , document , sum_cosine_map , all_cos_sum_queries );

	if(proba != 0){
		
		res += log( proba );

	}

	return res; 

}



//Same as before but over the entire collection
unordered_map <int,double> Dirichlet_embedding_model_Qexp(const double &mu , const vector<string> &query , const int query_id , const unordered_map< int ,  vector<string> > &collection , const unordered_map <string,int>  &cf , const unordered_map<string , double> &sum_cosine_map , const unordered_map< string , unordered_map<string,double> > &cosine_map , unordered_map< int , unordered_map<string,double> > &all_cos_sum_queries , const int collection_size ){
	
	unordered_map <int,double> list_doc;

	double proba;

	auto iterator = collection.begin();

	while(iterator != collection.end()){

		proba = Dirichlet_embedding_model_Qexp(mu , query , query_id , iterator->second , cf , sum_cosine_map , cosine_map , all_cos_sum_queries ,collection_size );
		if(proba!=0){list_doc[iterator->first]=proba;}

		iterator++;
		
	}

	return list_doc;

}



//Same as before but with all the queries and sort documents by their score
vector< vector< pair<int,double> > > Dirichlet_embedding_model_Qexp(const double &mu , const unordered_map< int ,  vector<string> > &queries , const unordered_map< int , vector<string> > &collection , const unordered_map <string,int>  &cf , const unordered_map<string , double> &sum_cosine_map , const unordered_map< string , unordered_map<string,double> > &cosine_map , unordered_map< int , unordered_map<string,double> > &all_cos_sum_queries , const int k , const int collection_size){

	vector< vector< pair<int,double> > > list_docs(queries.size());
	
	auto iterator = queries.begin();

	while(iterator != queries.end()){

		list_docs[iterator->first] =  kfirst_docs( Dirichlet_embedding_model_Qexp(mu , iterator->second , iterator->first , collection , cf , sum_cosine_map , cosine_map , all_cos_sum_queries , collection_size) , k );

		iterator++;

	}

	cout<<endl;

	return list_docs;

}






//Same as before but over the entire collection
unordered_map <int,double> FULL_embedding_model( const vector<string> &query , const unordered_map< int , vector<string> > &collection , Embedding &embedding ){
	
	unordered_map <int,double> list_doc;

	double proba;

	auto iterator = collection.begin();
	
	while(iterator != collection.end()){

		proba = embedding.cosine(query , iterator->second); 
		if(proba!=0){list_doc[iterator->first]=proba;}

		iterator++;

	}

	return list_doc;

}


//Same as before but with all the queries and sort documents by their score
vector< vector< pair<int,double> > > FULL_embedding_model( const unordered_map< int ,  vector<string> > &queries , const unordered_map< int ,  vector<string> > &collection , Embedding &embedding , const int k ){

	vector< vector< pair<int,double> > > list_docs(queries.size());

	auto iterator = queries.begin();
	
	while(iterator != queries.end() ){

		//cout<<"\rProcessing query "<<i+1<<"/"<<queries.size()<<flush;

		list_docs[iterator->first] = kfirst_docs( FULL_embedding_model( iterator->second , collection , embedding) , k );

		iterator++;

	}

	cout<<endl;

	return list_docs;

}



//Write the results in a file
void write_res_file(const vector< vector< pair<int,double> > > &results , const string &file_name , const string &query_Id , const double &lambda){

	ofstream myfile;
  	myfile.open (file_name.c_str() , ofstream::trunc);

	int rank;

	for(unsigned int i = 0 ; i < results.size() ; i++){
	
		rank = 0;
		if(i < 25){

			for(unsigned int j = 0 ; j < results[i].size() ; j++){

				myfile.precision(17);
				myfile << query_Id << setfill('0') << setw(3) << i+1 << " Q0 "<< results[i][j].first <<" "<< rank <<" "<< results[i][j].second <<" Hiemstra_LM"<<lambda<<"\n";
				rank++;

			}


		}

		else{

			for(unsigned int j = 0 ; j < results[i].size() ; j++){

				myfile.precision(17);
				myfile << query_Id << setfill('0') << setw(3) << i+26 << " Q0 "<< results[i][j].first <<" "<< rank <<" "<< results[i][j].second <<" Hiemstra_LM"<<lambda<<"\n";
				rank++;

			}


		}
	}

    myfile.close();

}

//Read a result file produced by trec_eval and return the map
vector<string> read_res_file(const string &file_name){

	FILE* f = fopen(file_name.c_str(),"r");

	size_t len = 10000;

	char* line = (char*)malloc(len);

	int current_line = 1;

	size_t nread;

	vector<string> res;
	string fra;
	cout<<"COUCOU"<<endl;
	int query_id = 1;
	while(int(nread = getline(&line, &len, f)) != -1 ){
		
		if(current_line%27==4){

			fra = string(line);
			fra = fra.substr(fra.size()-7,6);
			res.push_back( fra );
			query_id++;
		}
		current_line++;

	}

	free(line);

	return res;

}



//Write a <string,double> map in a file
void write_Average_Precision(const vector< string> &map , const string &file_name){

	ofstream myfile;

  	myfile.open(file_name.c_str() , ofstream::trunc);

	for(unsigned int i = 0 ; i < map.size() ; i++){

		myfile << map[i] + " , "  ;
	
	}

  	myfile.close();

}


/*
//Read a result file produced by trec_eval and return the map
vector<string> read_res_files(const string &directory_name){



	DIR *dir;
	struct dirent *ent;

	vector<string> set_files;

	if ((dir = opendir (directory_name.c_str())) != NULL) {


	  while ((ent = readdir (dir)) != NULL) {

		if(string(ent->d_name)!="." && string(ent->d_name)!=".."){set_files.push_back(directory_name + string(ent->d_name));}

	  }

	  closedir (dir);

	} 

	else {

	  cout<<"IMPOSSIBLE TO OPEN DIRECTORY"<<endl;

	}

	for(unsigned int i = 0 ; i < set_files.size() ; i++){

		set_files[i] = set_files[i] + " map " + read_res_file(set_files[i]);
	
	}

	for(unsigned int i = 0 ; i < set_files.size() ; i++){

		cout<<set_files[i]<<endl;
	
	}

	return set_files;
}
*/


//Read the collection and the set of queries and write the collection frequency and document frequency of the terms in files 
void write_all_files(const string &collection_file , const string &queries_file , const string &cf_file , const string &df_file){

	unordered_map< int ,  vector<string> > collection = read_file(collection_file);
	unordered_map< int ,  vector<string> > queries = read_file(queries_file);
	unordered_map <string,int>  cf = build_cf(collection);
	write_tf_file(cf , cf_file);
	unordered_map <string,int> df = build_df(collection, cf);
	write_tf_file(df , df_file);

}


//Builds the collection, queries cf and df
void get_all_info(const string &collection_file , const string &queries_file , unordered_map< int , vector<string> > &collection , unordered_map< int , vector<string> > &queries , unordered_map <string,int> &cf , unordered_map <string,int> &df){

	collection = read_file(collection_file);
	//collection = read_some_file(collection_file);
	queries = read_file(queries_file);
	cf = build_cf(collection);
	df = build_df(collection, cf);

}



//Read the similarities file without considering the similarities lower than a threshold
void delete_low_similarities(const string &collection_cosine_file , const string &queries_cosine_file , unordered_map< string , unordered_map<string,double> > &all_cos , unordered_map<string,double> &all_sum_cos , const double &threshold){

	unordered_map< string , unordered_map<string,double> > queries_cosine =  read_cos_map_map_file(queries_cosine_file , threshold);
	
	//unordered_map<string , double> all_sum_cos_queries  = read_cos_map_file(sum_queries_cosine_file , threshold);

	all_cos =  read_cos_map_map_file(collection_cosine_file , threshold);

	//all_sum_cos  = read_cos_map_file(sum_collection_cosine_file , threshold);

	aggregate_map( all_cos , queries_cosine );

	//aggregate_map( all_sum_cos , all_sum_cos_queries );

	all_sum_cos = all_fast_sum_cos_sim(all_cos);
	
}


void get_all_info2( const string &collection_file , const string &queries_file , const string &collection_cosine_file , const string &queries_cosine_file , unordered_map< int , vector<string> > &collection , unordered_map< int ,  vector<string> > &queries , unordered_map <string,int> &cf , unordered_map <string,int> &df , unordered_map< string , unordered_map<string,double> > &all_cos , unordered_map<string,double> &all_sum_cos , const double &threshold){

	collection = read_file(collection_file);
	queries = read_file(queries_file);
	cf = build_cf(collection);
	df = build_df(collection, cf);

	cout<<"Location of the cosine files :"<< queries_cosine_file <<" and "<< collection_cosine_file <<endl;

	clock_t begin = clock();
	
	unordered_map< string , unordered_map<string,double> > queries_cosine =  read_cos_map_map_file(queries_cosine_file , threshold);
	
	//unordered_map<string , double> all_sum_cos_queries  = read_cos_map_file(sum_queries_cosine_file , threshold);

	all_cos =  read_cos_map_map_file(collection_cosine_file , threshold);

	//all_sum_cos  = read_cos_map_file(sum_collection_cosine_file , threshold);

	aggregate_map( all_cos , queries_cosine );

	//aggregate_map( all_sum_cos , all_sum_cos_queries );

	all_sum_cos = all_fast_sum_cos_sim(all_cos);

	clock_t end = clock();
  	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<< "Time to read the cosine files : "<< elapsed_secs <<endl;
	

	
}


void get_all_info3( const string &collection_file , const string &queries_file , const string &collection_cosine_file , const string &queries_cosine_file , const string &sum_queries_cosine_file , unordered_map< int , vector<string> > &collection , unordered_map< int , vector<string> > &queries , unordered_map <string,int> &cf , unordered_map <string,int> &df , unordered_map< string , unordered_map<string,double> > &all_cos , unordered_map<string,double> &all_sum_cos , unordered_map< int , unordered_map<string,double> > &all_cos_sum_queries , const double &threshold){

	collection = read_file(collection_file);
	queries = read_file(queries_file);
	cf = build_cf(collection);
	df = build_df(collection, cf);

	cout<<"Location of the cosine files :"<< queries_cosine_file <<" and "<< collection_cosine_file <<endl;

	clock_t begin = clock();
	
	unordered_map< string , unordered_map<string,double> > queries_cosine =  read_cos_map_map_file(queries_cosine_file , threshold);
	
	//unordered_map<string , double> all_sum_cos_queries  = read_cos_map_file(sum_queries_cosine_file , threshold);

	all_cos =  read_cos_map_map_file(collection_cosine_file , threshold);

	//all_sum_cos  = read_cos_map_file(sum_collection_cosine_file , threshold);

	aggregate_map( all_cos , queries_cosine );

	//aggregate_map( all_sum_cos , all_sum_cos_queries );

	all_sum_cos = all_fast_sum_cos_sim(all_cos);

	all_cos_sum_queries =  read_cos_sum_query_map_map_file(sum_queries_cosine_file , threshold);

	clock_t end = clock();
  	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<< "Time to read the cosine files : "<< elapsed_secs <<endl;
	
}



//Reads everything
void fast_setup(const string &collection_file , const string &queries_file , const string &cf_file , const string &df_file , unordered_map< int , vector<string> > &collection , unordered_map< int , vector<string> > &queries , unordered_map <string,int> &cf , unordered_map <string,int> &df){

	collection = read_file(collection_file);
	queries = read_file(queries_file);
	cf = read_tf_file(cf_file);
	df = read_tf_file(df_file);

}


//Displays all the words in the collection that have a higher frequency than top
int display_top_frequent_words(const int top , const unordered_map <string,int>  &cf){

	testit p = cf.begin();

	int nb_words = 0;

	while( p!=  cf.end() ){
		
		nb_words += p->second;

		if(p->second > top){

			cout<<"Word of the vocabulary and it's frequency : ( "<< p->first <<" , "<< p->second <<" )"<<endl;

		}		
		p++;

	}

	cout<<"Sum of all the values : " << nb_words <<endl;

	return nb_words;

}



//display stuff
size_t display_stuff(const unordered_map< int , vector<string> > &collection , const unordered_map< int , vector<string> > &queries , const unordered_map <string,int> &cf , const unordered_map <string,int> &df){


	size_t queries_size = 0;

	auto iterator = queries.begin();

	while( iterator != queries.end() ){

		queries_size += iterator->second.size();

		iterator++;

	}
	cout<< "Size of the collection : "<< collection.size() <<endl;
	cout<<"Number of queries : "<< queries.size() <<endl;
	cout<<"Number of words in the queries : "<< queries_size <<endl;	
	cout<<"Size of the vocabulary : "<< cf.size() <<endl;
	cout<<"Size of the vocabulary : "<< df.size() <<endl;
	
	size_t nb_words = 0;
	testit p = cf.begin();

	while( p!=  cf.end() ){
		nb_words += p->second;
		p++;
	}

	cout<<"Number of words in the collection : " << nb_words <<endl;
	cout<<"**********************************************************"<<endl;
	cout<<"First document : "<<endl;
	display(collection.find(1103048)->second);
	cout<<"**********************************************************"<<endl;
	cout<<"First document : "<<endl;
	display(collection.find(1103049)->second);
	cout<<"**********************************************************"<<endl;
	cout<<"First document : "<<endl;
	display(collection.find(1103051)->second);
	cout<<"**********************************************************"<<endl;
	cout<<"First document : "<<endl;
	display(collection.find(1103054)->second);
	cout<<"**********************************************************"<<endl;
	cout<<"Neighborhood of size 3 at position 8 : "<<endl;
	display( neighbors(collection.find(0)->second , 8 , 3) );
	cout<<"**********************************************************"<<endl;
	cout<<"Last document : "<<endl;
	display(collection.find(collection.size()-1)->second);
	cout<<"**********************************************************"<<endl;

	return nb_words;

}

//return the size of the collection
size_t get_size_collection( const unordered_map <string,int> &cf ){

	size_t nb_words = 0;
	testit p = cf.begin();

	while( p!=  cf.end() ){
		nb_words += p->second;
		p++;
	}
	return nb_words;

}



//return the size of the collection
size_t get_size_flickr_collection( const unordered_map <string,int> &cf ){

	size_t nb_words = 0;

	testit p = cf.begin();

	while( p!=  cf.end() ){

		if(p->second != -1){
			
			nb_words += p->second;

		}

		p++;
	}

	return nb_words;

}




//count and displays the number of words in the vocabulary that have an embedding 
int nb_embedded_words_in_voc(Embedding &embedding , const unordered_map <string,int> &cf){

	int count = 0;
	int nb_words_tot = 0;
	int	nb_words_emb = 0;
	int nb_iter = 0;

	testit p = cf.begin();

	while( p!=  cf.end() ){
		
		if( embedding.get(p->first.c_str()) != nullptr ){

			count++;
			nb_words_emb += p->second;

		}

		nb_words_tot += p->second;		
		nb_iter++;
		p++;

	}
	
	cout<<"Percentage of words of the vocabulary that have an embedding : "<< double(count)*100/cf.size()<<"%"<<endl;
	cout<<"Percentage of words of the collection that have an embedding : "<< double(nb_words_emb)*100/nb_words_tot<<"%"<<endl;

	return count;

}



//count the number of words in the queries that have an embedding 
int nb_embedded_words_in_queries(Embedding &embedding , const unordered_map< int , vector<string> > &queries , const bool displ){

	unsigned int count = 0;
	unsigned int nb_words_tot = 0;
	unsigned int temp = 0;
	unsigned int nb_empty_queries = 0;
	unsigned int nb_full_queries = 0;

	auto iterator = queries.begin();

	while(iterator != queries.end()){

		temp = iterator->second.size();

		for(unsigned int j = 0 ; j < iterator->second.size() ; j++){

			if( embedding.get(iterator->second[j].c_str()) != nullptr ){
				count++;
				temp--;

			}

			else if(displ==true){

					//cout<<queries[i][j]<<endl;

				}

			nb_words_tot++;
			
		}
		
		if(temp == iterator->second.size()){

			nb_empty_queries++;
			if(displ == true){
				cout<<"Empty Query "<<iterator->first + 1<<" ";
				display(iterator->second);
			}

		}


		if(temp == 0){

			nb_full_queries++;
			if(displ == true){
				cout<<"Full Query "<<iterator->first + 1<<" ";
				display(iterator->second);
			}

		}
	
		iterator++;

	}
	
	cout<<"Percentage of words of the queries that have an embedding : "<< double(count)*100/nb_words_tot<<"%"<<endl;
	cout<<"Percentage of queries that does not have any term that have an embedding : "<< double(nb_empty_queries)*100/queries.size()<<"%"<<endl;
	cout<<"Percentage of queries that have all their terms that have an embedding : "<< double(nb_full_queries)*100/queries.size()<<"%"<<endl;

	return count;

}




//Performs a set of experiments with a "regular" model
void launch_Hiemstra_experience(const string &collection_file , const string &queries_file , double lambda , const double lambda_step , const int k){

	vector< vector< pair<int,double> > > results;
	unordered_map< int , vector<string> > collection;
	unordered_map< int , vector<string> > queries;
	unordered_map <string,int> cf;
	unordered_map <string,int> df;
	string file_name;

	clock_t begin = clock();

	get_all_info(collection_file , queries_file , collection , queries , cf , df);

	write_xml_collection_file(collection, "../../../Terrier4.2/terrier-core-4.2/share/CHIC2012/synthetic/xml_synthetic_collection");

	display_stuff(collection ,queries , cf , df);

	clock_t end = clock();
  	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<< "Time to do all the stuff : "<< elapsed_secs <<endl;

	display_stuff(collection , queries , cf , df);

	size_t nb_words = get_size_collection(cf);

	while(lambda < 1 + lambda_step/1000){

		begin = clock();
		results = Hiemstra_language_model(queries , collection , cf , nb_words , k , lambda);
		end = clock();
		elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		cout<< "Time to perform all the queries for lambda = "<<lambda<<" : "<< elapsed_secs <<endl;
		file_name = "../data/res/logresults";
		file_name += to_string(lambda); 
		write_res_file(results , file_name , "CHIC-" , lambda);
		lambda += lambda_step;

	}

}



//Performs a set of experiments with an embedded model
void launch_Dirichlet_experience(const string &collection_file , const string &queries_file , double &mu , const double &mu_step , const int nb_iter , const int k ){

	vector< vector< pair<int,double> > > results;
	unordered_map< int , vector<string> > collection;
	unordered_map< int , vector<string> > queries;
	unordered_map <string,int> cf;
	unordered_map <string,int> df;
	string file_name;
	double mu_temp;
	clock_t begin = clock();

	get_all_info( collection_file , queries_file , collection , queries , cf , df);

	write_xml_collection_file(collection, "../../../Terrier4.2/terrier-core-4.2/share/CHIC2012/synthetic/xml_synthetic_collection");

	clock_t end = clock();
  	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<< "Time to do all the stuff : "<< elapsed_secs <<endl;

	display_stuff(collection , queries , cf , df);

	size_t nb_words = get_size_collection(cf);

	int nthreads, tid;

	#pragma omp parallel private( file_name , mu_temp , results , tid) shared(mu , mu_step) 
	{

		tid = omp_get_thread_num();

        //cout<<"Thread No "<<tid<<endl;
        if (tid == 0)
        {
            nthreads = omp_get_num_threads();
            cout<<"Number of thread: "<<nthreads<<endl;
        }


		#pragma omp for schedule(static)
		for(int current_iter = 0 ; current_iter < nb_iter ; current_iter++){

			mu_temp = mu + mu_step*(current_iter);
			results = Dirichlet_language_model(mu_temp , queries , collection , cf , k , nb_words);
			file_name = "../data/res/results";
			file_name += to_string(mu_temp);
			write_res_file(results , file_name , "CHIC-" , mu_temp);
			//mu += mu_step;
			//current_iter++;
		
		}

	}

}




//Performs a set of experiments with an embedded model
void launch_embedded_experience(const string &collection_file , const string &queries_file , const string &collection_cosine_file , const string &queries_cosine_file , double &mu , const double &mu_step , const int nb_iter , const int k , double &threshold , const double &threshold_min , const double &threshold_max , const double &threshold_step , const double &alpha , const double &alpha_min , const double &alpha_max , const double &alpha_step){

	unordered_map< string , unordered_map<string,double> > all_cos;
	unordered_map<string,double> all_sum_cos;
	vector< vector< pair<int,double> > > results;
	unordered_map< int , vector<string> > collection;
	unordered_map< int , vector<string> > queries;
	unordered_map <string,int> cf;
	unordered_map <string,int> df;
	string file_name;
	double mu_temp;
	double alpha_temp;
	clock_t begin = clock();

	get_all_info2( collection_file , queries_file , collection_cosine_file , queries_cosine_file , collection , queries , cf , df , all_cos , all_sum_cos , threshold);

	write_xml_collection_file(collection, "../../../Terrier4.2/terrier-core-4.2/share/CHIC2012/synthetic/xml_synthetic_collection");

	clock_t end = clock();
  	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<< "Time to do all the stuff : "<< elapsed_secs <<endl;

	display_stuff(collection , queries , cf , df);

	size_t nb_words = get_size_collection(cf);

	//int current_iter;

	int nthreads, tid;
		
	while(threshold < threshold_max + abs(threshold_step/1000) && threshold > threshold_min - abs(threshold_step/1000)){

		//current_iter = 0;
		
		alpha_temp = alpha;

		while(alpha_temp < alpha_max + abs(alpha_step/1000) && alpha_temp > alpha_min - abs(alpha_step/1000)){

			#pragma omp parallel private(file_name , mu_temp , results , tid) shared(mu , mu_step , alpha_temp) 
			{

				tid = omp_get_thread_num();

				//cout<<"Thread No "<<tid<<endl;
				if (tid == 0)
				{
				    nthreads = omp_get_num_threads();
				    cout<<"Number of thread: "<<nthreads<<endl;
				}

				#pragma omp for schedule(static)
				for(int current_iter = 0 ; current_iter < nb_iter ; current_iter++){
				//while(current_iter < nb_iter){
					mu_temp = mu + mu_step*(current_iter);
					results = Dirichlet_embedding_model(mu_temp , queries , collection , cf , all_sum_cos , all_cos , k , nb_words , alpha_temp);
					cout<< "Performed all the queries for mu = "<<mu_temp<< " , for alpha = "<< alpha_temp <<" and the threshold = " << threshold <<endl;
					file_name = "../data/res/results";
					file_name += to_string(mu_temp);
					file_name += "thresh";
					file_name += to_string(threshold);
					file_name += "alpha";
					file_name += to_string(alpha_temp);			
					write_res_file(results , file_name , "CHIC-" , mu_temp);
					//mu_temp += mu_step;
					//current_iter++;

				}

			}

			alpha_temp += alpha_step;
	
		}

	threshold += threshold_step;

	delete_low_similarities(collection_cosine_file , queries_cosine_file , all_cos , all_sum_cos , threshold);

	}

}


//Performs a set of experiments with an embedded model
void launch_expanded_experience(const string &collection_file , const string &queries_file , const string &collection_cosine_file , const string &queries_cosine_file , const string &sum_queries_cosine_file , double &mu ,  const double &mu_step , const int nb_iter , const int k , double &threshold , const double &threshold_min , const double &threshold_max , const double &threshold_step ){

	unordered_map< string , unordered_map<string,double> > all_cos;
	unordered_map< int , unordered_map<string,double> > all_cos_sum_queries;
	unordered_map<string,double> all_sum_cos;
	vector< vector< pair<int,double> > > results;
	unordered_map< int , vector<string> > collection;
	unordered_map< int , vector<string> > queries;
	unordered_map <string,int> cf;
	unordered_map <string,int> df;
	string file_name;
	double mu_temp;


	clock_t begin = clock();

	get_all_info3( collection_file , queries_file , collection_cosine_file , queries_cosine_file , sum_queries_cosine_file , collection , queries , cf , df , all_cos , all_sum_cos , all_cos_sum_queries , threshold);

	write_xml_collection_file(collection, "../../../Terrier4.2/terrier-core-4.2/share/CHIC2012/synthetic/xml_synthetic_collection");

	clock_t end = clock();
  	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<< "Time to do all the stuff : "<< elapsed_secs <<endl;

	display_stuff(collection , queries , cf , df);

	size_t nb_words = get_size_collection(cf);

	//int current_iter;

	int nthreads, tid;
		
	while(threshold < threshold_max + abs(threshold_step/1000) && threshold > threshold_min - abs(threshold_step/1000)){

		//current_iter = 0;

	

	#pragma omp parallel private(file_name , mu_temp , results , tid) shared(mu , mu_step ) 
	{
		

		tid = omp_get_thread_num();

        //cout<<"Thread No "<<tid<<endl;
        if (tid == 0)
        {
            nthreads = omp_get_num_threads();
            cout<<"Number of thread: "<<nthreads<<endl;
        }

		mu_temp = mu;

		#pragma omp for schedule(static)
		for(int current_iter = 0 ; current_iter < nb_iter ; current_iter++){
		//while(current_iter < nb_iter){
			mu_temp = mu + mu_step*(current_iter);
			results = Dirichlet_embedding_model_Qexp(mu_temp , queries , collection , cf , all_sum_cos , all_cos , all_cos_sum_queries , k , nb_words);
			cout<< "Performed all the queries for mu = "<<mu_temp<< " and the threshold = " <<threshold<<endl;
			file_name = "../data/res/results";
			file_name += to_string(mu_temp);
			file_name += "thresh";
			file_name += to_string(threshold);			
			write_res_file(results , file_name , "CHIC-" , mu_temp);
			//mu_temp += mu_step;
			//current_iter++;

		}

	}

	threshold += threshold_step;

	delete_low_similarities(collection_cosine_file , queries_cosine_file , all_cos , all_sum_cos , threshold);

	}


}



//Performs a set of experiments with an embedded model
void launch_full_embedding_experience(const string &collection_file , const string &queries_file , const string &embeddings_file , const int k ){

	vector< vector< pair<int,double> > > results;
	unordered_map< int , vector<string> > collection;
	unordered_map< int , vector<string> > queries;
	unordered_map <string,int> cf;
	unordered_map <string,int> df;
	Embedding embedding;

	embedding.load_Word2VecBinFormat(embeddings_file.c_str());

	clock_t begin = clock();

	get_all_info( collection_file , queries_file , collection , queries , cf , df);


	display_stuff(collection , queries , cf , df);

	clock_t end = clock();
  	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<< "Time to read the files : "<< elapsed_secs <<endl;

	begin = clock();

	//results = FULL_embedding_model( queries , collection , embedding , k );
	FULL_embedding_model( queries[0] , collection , embedding);

	end = clock();
  	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<< "Time to compute the results : "<< elapsed_secs <<endl;

	//write_res_file(results , "../data/res/results" , "CHIC-" , 0);


}























void check_flickr_string(string &s){

	for(unsigned int i = 0 ; i < s.size() ; i++){

		if( !isValidChar(s[i]) ){
					
			s.erase(s.begin()+i);
			i--;					

		}

	}

}


//Reads a collection/set of queries in the input file
unordered_map< string , vector<string> > read_queries_flickr_file(const string &file_name){

 	FILE* f = fopen(file_name.c_str(),"r");

	size_t len = 10000;

	char* line = (char*)malloc(len);

	size_t nread;
	
	vector<string> temp;

	string temp2;

	unordered_map< string , vector<string> > queries;

	while (int(nread = getline(&line, &len, f)) != -1) {

		temp = easy_split(string(line) , ' ');

		if(temp.size() > 2){

			temp2 = temp[0];

			temp.erase(temp.begin() , temp.begin() + 2);

			for(unsigned int i = 0 ; i < temp.size() ; i++){

				if(temp[i] == "\n" || temp[i].size() == 0 ){
					
					temp.erase(temp.begin()+i);
					i--;					

				}

			}
	
			queries[temp2] = temp;

		}
		
	}

	free(line);	

	return queries;
	
}



//Reads a collection/set of queries in the input file
unordered_map< string , vector<string> > read_collection_flickr_file(const string &file_name){

 	FILE* f = fopen(file_name.c_str(),"r");

	size_t len = 10000;

	char* line = (char*)malloc(len);

	size_t nread;
	
	vector<string> temp;

	vector<string> previous_temp;

	string temp2;

	unordered_map< string , vector<string> > collection;

	while (int(nread = getline(&line, &len, f)) != -1) {

		//cout<<"Begin line"<<line<<"end line"<<endl;

		temp = easy_split(string(line) , ' ');

		//cout<<temp.size()<<endl;
		
		if(temp.size() > 2){

			//cout<<"Inside if"<<endl;

			
			if(temp[0].find("/video/webvision/flickr/") != string::npos){

				temp2 = temp[0];

				temp.erase(temp.begin() , temp.begin() + 2);

				for(unsigned int i = 0 ; i < temp.size() ; i++){

					if(temp[i] == "\n" || temp[i].size() == 0 ){
					
						temp.erase(temp.begin()+i);
						i--;					

					}

				}


				for(unsigned int i = 0 ; i < temp.size() ; i++){

					check_flickr_string(temp[i]);

				}

		
				collection[temp2] = temp;

				previous_temp = temp;

			}
			
			else{

				temp.erase(temp.begin() , temp.begin() + 2);

				for(unsigned int i = 0 ; i < temp.size() ; i++){

					if(temp[i] == "\n" || temp[i].size() == 0 ){
					
						temp.erase(temp.begin()+i);
						i--;					

					}

				}


				for(unsigned int i = 0 ; i < temp.size() ; i++){

					check_flickr_string(temp[i]);

				}			


				temp.insert( temp.end(), previous_temp.begin(), previous_temp.end() );


				collection[temp2] = temp;

				previous_temp = temp;
		
			}

		}
		
	}

	free(line);	

	return collection;
	
}






//Builds the vocabulary given the collection in input and stores it in an unordered_map of <string,int> the string being the words and the int their associated collection frequency
unordered_map <string,int> build_voc_flickr(const unordered_map< string , vector<string> > &collection , const unordered_map< string , vector<string> > &queries){

	unordered_map <string,int>  cf;

	auto iterator = collection.begin();

	while( iterator != collection.end() ){

		for(unsigned int j = 0 ; j < iterator->second.size() ; j++){
					
			if(iterator->second[j].size() > 1 || isValidChar(iterator->second[j][0])){
				cf[ iterator->second[j] ]++;
			}
		}

		iterator++;
	
	}

	iterator = queries.begin();

	while( iterator != queries.end() ){

		for(unsigned int j = 0 ; j < iterator->second.size() ; j++){
					
			if(iterator->second[j].size() > 1 || isValidChar(iterator->second[j][0])){

				if(cf.find(iterator->second[j]) == cf.end()){cf[ iterator->second[j] ] = -1;}
			
			}
		}

		iterator++;
	
	}


	return cf;
	
}


//Same as before but over the queries
unordered_map< string , unordered_map<string,double> > closest_terms(const unordered_map< string , vector<string> > &queries , const unordered_map <string,int> &cf , Embedding &embedding ,  const double &threshold){

	unordered_map< string , unordered_map<string,double> > set_most_sim;

	int compteur = 0;

	auto iterator = queries.begin();

	while(iterator != queries.end()){

		for(unsigned int j = 0 ; j < iterator->second.size() ; j++){

			set_most_sim[iterator->second[j]] = closest_terms(iterator->second[j] , cf , embedding , threshold);

		}

		compteur++;

		cout<<"\rProgress : ("<< compteur <<"/" << queries.size() << ")"<<flush;

		iterator++;

	}

	return set_most_sim;

}



//Computes the log of the probability that query was generated by the document w.r.t dirichlet language model  pdir(query | document)
double flickr_embedding_model(const double &mu , const vector<string> &query , const vector<string> &document , const unordered_map <string,int>  &cf , const unordered_map<string , double> &sum_cosine_map , const unordered_map< string , unordered_map<string,double> > &cosine_map , const int collection_size , const double &alpha ){

	double res = 0;

	double proba;

	int collection_freq;

	for(unsigned int i = 0 ; i < query.size() ; i++){

		
		proba = proba_doc_generate_term(query[i] , document , sum_cosine_map , cosine_map , alpha);

		collection_freq = coll_freq(cf , query[i]);

		if(collection_freq != 0 && proba != 0){

			res += log( ( (proba*document.size()) + mu*((double)collection_freq/collection_size))/(document.size() + mu) );
		
		}

		else if(collection_freq != 0 && proba == 0){

			res += log( (double)collection_freq/collection_size );

		}

		else if(collection_freq == 0 && proba != 0){

			res += log( proba );

		}

	}

	return res; 

}



int main(int argc, char** argv) {

	cout<<"Program start"<<endl;

	if(argc == 1){

		//Read a result file produced by trec_eval and return the map
		vector<string> Average_Precision = read_res_file("../../../trec_eval.8.1/results/resultat0");

		cout<<"Taille vecteur : "<< Average_Precision.size() <<endl;

 		write_Average_Precision(Average_Precision , "../data/AveragePrecision/AP");

	}

	string collection_file = "../data/collection/porter_stop_string_content";
	string queries_file = "../data/queries/porter_stop_queries2012.txt";

	if(argc > 1 && string(argv[1]) == "cosine"){

		string collection_cosine_file = "../data/embeddings/GoogleNews-vectors-negative300/cos04";
		string queries_cosine_file = "../data/embeddings/GoogleNews-vectors-negative300/cos_queries_04";
		string sum_queries_cosine_file = "../data/embeddings/GoogleNews-vectors-negative300/cos_sum_queries_04";
		string embeddings_file = "../data/embeddings/GoogleNews-vectors-negative300/GoogleNews-vectors-negative300.bin";
		//string sum_proba_file = "../data/embeddings/GoogleNews-vectors-negative300/sum_proba04";


		if(argc > 2 && string(argv[2]) == "google"){

			collection_cosine_file = "../data/embeddings/GoogleNews-vectors-negative300/cos04";
			queries_cosine_file = "../data/embeddings/GoogleNews-vectors-negative300/cos_queries_04";
			sum_queries_cosine_file = "../data/embeddings/GoogleNews-vectors-negative300/cos_sum_queries_04";
			embeddings_file = "../data/embeddings/GoogleNews-vectors-negative300/GoogleNews-vectors-negative300.bin";
			//sum_proba_file = "../data/embeddings/GoogleNews-vectors-negative300/sum_proba04";

		}

		if(argc > 2 && string(argv[2]) == "wiki"){

			collection_cosine_file = "../data/embeddings/wiki/cos04";
			queries_cosine_file = "../data/embeddings/wiki/cos_queries_04";
			sum_queries_cosine_file = "../data/embeddings/wiki/cos_sum_queries_04";
			embeddings_file = "../data/embeddings/wiki/wiki.bin";
			//sum_proba_file = "../data/embeddings/wiki/sum_proba04";

		}

		Embedding embedding;

		embedding.load_Word2VecBinFormat(embeddings_file.c_str());

		unordered_map< int , vector<string> > collection;

		unordered_map< int , vector<string> > queries;

		unordered_map <string,int> cf;

		unordered_map <string,int> df;

		get_all_info( collection_file , queries_file , collection , queries , cf , df);

		cout<<"Query size : "<< queries.size() <<endl;

		nb_embedded_words_in_voc(embedding , cf);

		nb_embedded_words_in_queries(embedding , queries , false);

		unordered_map< string , unordered_map<string,double> > set_closest_words;

		unordered_map< int , unordered_map<string,double> > set_closest_words_sum_query;

		return 0;

		//Queries
		
		set_closest_words = closest_terms(queries , cf ,embedding ,  0.4);

		write_map_map(set_closest_words , queries_cosine_file);

		//Sum Queries

		set_closest_words_sum_query = closest_terms_sum_query(queries , cf , embedding , 0.4);

		write_map_map(set_closest_words_sum_query , sum_queries_cosine_file);

		//Vocabulary

		set_closest_words = closest_terms(cf , embedding , 0.4);

		write_map_map(set_closest_words , collection_cosine_file);

		

	}

	else if(argc > 1 && string(argv[1]) == "hiemstra"){


		double lambda = 0.0;
		int k = 1000;
		double lambda_step = 0.05;
		launch_Hiemstra_experience(collection_file , queries_file , lambda , lambda_step , k);

	}

	else if(argc > 1 && string(argv[1]) == "dirichlet"){


		double mu = 0;
		int k = 1000;
		double mu_step = 20;
		int nb_iter = 1;
		launch_Dirichlet_experience(collection_file , queries_file , mu , mu_step , nb_iter , k );

	}


	else if(argc > 1 && string(argv[1]) == "embedding"){

		string collection_cosine_file = "../data/embeddings/GoogleNews-vectors-negative300/cos04";
		string queries_cosine_file = "../data/embeddings/GoogleNews-vectors-negative300/cos_queries_04";
		//string sum_proba_file = "../data/embeddings/GoogleNews-vectors-negative300/sum_proba04";
		

		if(argc > 2 && string(argv[2]) == "google"){

			collection_cosine_file = "../data/embeddings/GoogleNews-vectors-negative300/cos04";
			queries_cosine_file = "../data/embeddings/GoogleNews-vectors-negative300/cos_queries_04";
			//sum_proba_file = "../data/embeddings/GoogleNews-vectors-negative300/sum_proba04";

		}

		if(argc > 2 && string(argv[2]) == "wiki"){

			collection_cosine_file = "../data/embeddings/wiki/cos04";
			queries_cosine_file = "../data/embeddings/wiki/cos_queries_04";
			//sum_proba_file = "../data/embeddings/wiki/sum_proba04";

		}

		double mu = 130;

		int k = 1000;

		double mu_step = 10;

		int nb_iter = 1;

		double threshold = 0.7;

		double threshold_step = 1.0;

		double threshold_min = 0.0;

		double threshold_max = 1.0;

		double alpha = 0.0;
			
		double alpha_step = 0.1;

		double alpha_min = 0.0;

		double alpha_max = 1;

		launch_embedded_experience(collection_file , queries_file , collection_cosine_file , queries_cosine_file , mu , mu_step , nb_iter , k , threshold , threshold_min , threshold_max , threshold_step , alpha , alpha_min , alpha_max, alpha_step);

	}

	return 0;

}
