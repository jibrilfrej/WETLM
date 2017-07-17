//  Copyright 2013 Google Inc. All Rights Reserved.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h> // mac os x
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cassert>
#include <unordered_map>




const long long max_size = 2000;         // max length of strings
const long long N = 1;                   // number of closest words that will be shown
const long long max_w = 128;              // max length of vocabulary entries




using namespace std; 

typedef unordered_map<string,int>::const_iterator testit;

#define err1 0x0000 // la fonction lecture_fichier_bin n a pas trouve le repertoire en entree


//Check if the input char is a majuscule letter
inline
bool isMajChar(const char c){

	return ( c > 64 && c < 91 );

}


inline
bool isMinChar(const char c){

	return ( c > 96 && c < 123 );

}




 // A simple word
class Word {
private:
  char data[max_w];
public:
  Word(){}
  char& operator[](size_t i){return data[i];}
  bool operator!=(const char* w){return strcmp(data,w)!=0;}
  bool operator==(const char* w){return strcmp(data,w)==0;}
  const char* c_str()const{return data;}
  const string to_string()const{return string(data);}

};

inline
bool is_letter(char c){return (c >= 'a' && c <= 'z')||( c >= 'A' && c <= 'Z')||(c == '_')||(c >= '0' && c <= '9');}

class Sentence {
private:
  char *data;
  vector<char*> words;
  size_t max_size;

public:
  Sentence(size_t s=100*max_w):data((char*)malloc(s)),max_size(s){}
  ~Sentence(){free(data);}
  size_t size(){return words.size();}

  char* operator[](size_t i){

  assert(i < words.size());

  return words[i];
  }

  
  void read_file(FILE* f = stdin){
	unsigned int a = 0;

	while (1) {
      data[a] = fgetc(f);
      if ((data[a] == '\n') || (a >= max_size - 1)) {
        data[a] = 0;
        break;
      }
      a++;
    }
	
  }




void read_string(string s){
	unsigned int a = 0;

	while (1) {
      data[a] = s[a];
      if ((data[a] == '\n') || (a >= max_size - 1)) {
        data[a] = 0;
        break;
      }
      a++;
    }
	
  }



  void split(){
	
	char* t = data;
	while(*t != 0 && !is_letter(*t)){t++;}

	while(*t != 0){

		words.push_back(t);
		cout<<"boucle2 split : "<<t<<endl;
		while(*t != 0 && is_letter(*t)){t++;}
		if(*t != 0){
			*t = 0;
			t++;
		}		
		
		while(*t != 0 && !is_letter(*t)){t++;}

	}

	

  }
  
};

// To store a vocabulary and the embedded vectors
class Embedding {

public:

	Embedding(){}

	// Load all data from binary word2vect format
	// Return 0 if OK
	int load_Word2VecBinFormat(const char* filename ///< the name of the binary file
	);
    int load_Word2VecBinFormat(const std::string& filename ///< the name of the binary file
	) {return load_Word2VecBinFormat(filename.c_str());}
	
	//Save the vocabulary in a file
	void save_voc(const string &file_name);

	//Read the vocabulary from a file
	void set_hmap(const unordered_map<string,int> &stemmed , const string &file_name);

    // Return a ponter to the embedded vector of the word word
	float* get(const char* word);

	//return the number of terms in the query that has an embedding
	int check_embedding(const vector<string> &query );

	// Return the cosine similarity between two words
	float cosine(const char* , const char*);

	// Return the cosine similarity between a sums of words and a word
	float cosine(const vector<string>  &set1, const  string  &set2);

	// Return the cosine similarity between two sums of words
	float cosine(const vector<string>  &set1, const  vector<string>  &set2);

	// Write the expanded queries in the file named name_doc_expanded_queries
	void expand_queries(string name_doc_queries, string name_doc_expanded_queries);

	//Displays vocab and M sizes and the size of the embedded vectors 
	void display_attributes();

	//Get size of vocabulary
	size_t size_voc()const{return vocab.size();} 

	//Return the list of all the words (to use for debug only)
	const char* operator[](const unsigned int id ) const{return vocab[id].c_str();}


private:
   
	//hMAP of the vocabulary
	unordered_map<string,int> hmap; 

	// the table of all words
	vector<Word> vocab;

	// Matrix of all embeded vectors
	vector<float> M;

	//Size of each embedded vectors
	uint64_t size;

	// Return a ponter to the embedded vector that starts at indice i
	float* get(size_t i);

};


void Embedding::save_voc(const string &file_name){

	ofstream myfile;

  	myfile.open (file_name.c_str());

	auto iterator = hmap.begin();

	while(iterator != hmap.end()){

		myfile << iterator->first + " " + to_string(iterator->second) + "\n";

		iterator++;

	}

  	myfile.close();

}



//Take a string as input and return a pair of string and double using the delimiter
pair<string,int> read_line(const string &str, const char delimiter) {

  pair<string,int> res;
  stringstream ss(str); // Turn the string into a stream.
  string tok;
  tok[0] = 0;
  getline(ss, tok, delimiter);
  res.first = tok;
  getline(ss, tok, delimiter);
  res.second = atoi(tok.c_str());		
  
  return res;
}



//Read a file containing the map
void Embedding::set_hmap(const unordered_map<string,int> &stemmed , const string &file_name){

	FILE* f = fopen(file_name.c_str(),"r");

	size_t len = 10000;

	char* line = (char*)malloc(len);

	size_t nread;

	pair<string,int> temp;

	while (int(nread = getline(&line, &len, f)) != -1) {

		temp = read_line( string(line) , ' ');

		hmap[temp.first] = temp.second;
		
	}

	free(line);	

}








// Read in memory a binary file of vod2vec output
int Embedding::load_Word2VecBinFormat (
const char* file_name ///< the name of the binary file
	)
{
	int k = 0;
	cout<<"Before opening binary file"<<endl;
	// Open the binary file
	FILE *f = fopen(file_name, "rb");
	if(f == NULL) return int(err1);
	
	cout<<"After opening binary file"<<endl;

	uint64_t words;

	//Reading the header of the file 
	//Warning(1) : the integers must have a lower value than 2^64
	//Warning(2) : the scanf must get than l means 64 bits
	fscanf(f, "%lu", &words);
  	fscanf(f, "%lu", &size);
	
	cout<< "Size of the vocab : " <<words<<endl;
	cout<< "Size of the vectors : " <<size<<endl;
	// Reserve memory for the vocabulary 
	vocab.resize(words); 
	cout<<"After allocation of the vocab"<<endl;
	// Reserve memory for the matrix of embedded vectors 
	M.resize(words*size); 
	cout<<"After allocation of the vectors"<<endl;
	//separator char between the word and the embedded vector
	char ch;
	cout<<"Starting reading binary file : "<<endl;
	for (size_t b = 0; b < words; b++) 
	{
		
		cout<<"\rProcessing ...................."<< 100*(b+1)/words <<"%"<<flush;
		
		//Reads the bth word
		fscanf(f, "%s%c", &vocab[b][0], &ch);
		hmap[vocab[b].to_string()] = b;
		//Read the bth vector
		fread(&M[b * size], sizeof(float), size , f);
		//Normalization of the vectors
		float len = 0;
		for (size_t a = 0; a < size; a++) len += M[a + b * size] * M[a + b * size];
		len = sqrt(len);
		for (size_t a = 0; a < size; a++) M[a + b * size] /= len;
		for(size_t a = 0; a < size; a++){

			cout << M[a + b * size] <<endl;
			cin >> k;		
		} 
  	}
	cout<<endl;
	fclose(f);
	
	save_voc("non_stemmed_voc");

	return 0;
}







float* Embedding::get(const char* word){

	string term = string(word);

	auto it = hmap.find( term );

	if(it != hmap.end()){

		return get(it->second);

	}
	/*
	for( unsigned int i = 0 ; i < term.size() ; i++){

		if( term[i] > 47 && term[i] < 58 ){
		
			term[i]='#';
			
		}
	}

	it = hmap.find( term );

	if(it != hmap.end()){

		return get(it->second);

	}

	*/
	
	if( ( (term[0]) > 96 && (term[0]) < 123 ) ){

		term[0] = term[0] - 32;

		it = hmap.find( term );

		if(it != hmap.end()){

			return get(it->second);

		}

	}
		
	
	
	return nullptr;

}

inline
float* Embedding::get(size_t i){
	assert( i < vocab.size() );
	return &M[i*size];

}


int Embedding::check_embedding(const vector<string> &query){

	int res = 0;

	for(unsigned int i = 0 ; i < query.size() ; i++){

		if(get(query[i].c_str()) != nullptr){res++;}

	}

	return res;

}



float Embedding::cosine(const char* word1, const char* word2){

	float *vect1 = get(word1);
	float *vect2 = get(word2);
	if(vect1 == nullptr ||vect2 == nullptr){return 0;}
	float dist = 0;
	for(unsigned int i = 0; i < size ; i++){
	
		dist += (*vect1)*(*vect2);
		vect1++;
		vect2++;	

	}

	return dist;
}


float Embedding::cosine(const vector<string> &set1 , const string &set2 ){

	if(set1.size() == 0){return 0;}

	float *term2 = get(set2.c_str());

	if(term2==nullptr){return 0;}

	float *sumvect1 = nullptr;

	unsigned int iterator1 = 0;

	while(sumvect1 == nullptr && iterator1 < set1.size()){

		
		sumvect1 = get(set1[iterator1].c_str());
		iterator1++;

	}

	if(sumvect1 == nullptr){return 0;}

	float *vect;

	float len = 0;

	for(unsigned int i = iterator1+1 ; i < set1.size() ; i++){

		vect = get(set1[i].c_str());

		if(vect != nullptr){

				
				for(unsigned int j = 0; j < size ; j++){
	
					sumvect1[j] += (*vect);
					vect++;	

				}

		}
		
	}
	len = 0;
	for (unsigned int i = 0; i < size ; i++) len += sumvect1[i]*sumvect1[i];
	len = sqrt(len);
	for (unsigned int i = 0; i < size ; i++) sumvect1[i]/= len;

	
	float dist = 0;

	for(unsigned int i = 0; i < size ; i++){
	
		dist += (*sumvect1)*(*term2);
		sumvect1++;
		term2++;	

	}

	return dist;
	
}






float Embedding::cosine(const vector<string> &set1 , const vector<string> &set2 ){

	if(set1.size() == 0 || set2.size() == 0 ){return 0;}

	float *sumvect1 = nullptr;
	float *sumvect2 = nullptr;

	unsigned int iterator1 = 0;

	while(sumvect1 == nullptr && iterator1 < set1.size()){

		
		sumvect1 = get(set1[iterator1].c_str());
		iterator1++;

	}

	if(sumvect1 == nullptr){return 0;}

	unsigned int iterator2 = 0;

	while(sumvect2 == nullptr && iterator2 < set2.size()){

		
		sumvect2 = get(set2[iterator2].c_str());
		iterator2++;

	}

	if(sumvect2 == nullptr){return 0;}


	float *vect;

	float len = 0;

	for(unsigned int i = iterator1+1 ; i < set1.size() ; i++){

		vect = get(set1[i].c_str());

		if(vect != nullptr){

				
				for(unsigned int j = 0; j < size ; j++){
	
					sumvect1[j] += (*vect);
					vect++;	

				}

		}
		
	}
	len = 0;
	for (unsigned int i = 0; i < size ; i++) len += sumvect1[i]*sumvect1[i];
	len = sqrt(len);
	for (unsigned int i = 0; i < size ; i++) sumvect1[i]/= len;

	for(unsigned int i = 1 ; i < set2.size() ; i++){

		vect = get(set2[i].c_str());

		if(vect != nullptr){

			for(unsigned int j = 0; j < size ; j++){

				sumvect2[j] += (*vect);
				vect++;	

			}
				
		}
		
	}
	

	len = 0;
	for (unsigned int i = 0; i < size ; i++) len += sumvect2[i]*sumvect2[i];
	len = sqrt(len);
	for (unsigned int i = 0; i < size ; i++) sumvect2[i]/= len;

	
	float dist = 0;

	for(unsigned int i = 0; i < size ; i++){
	
		dist += (*sumvect1)*(*sumvect2);
		sumvect1++;
		sumvect2++;	

	}

	return dist;
	
}



//Displays vocab and M sizes and the size of the embedded vectors 
void Embedding::display_attributes(){

	cout<<"Size of vocab : "<<vocab.size()<<endl;
	cout<<"Size of M : "<<M.size()<<endl;
	cout<<"Size of vectors : "<<size<<endl;

}


void Embedding::expand_queries(string name_doc_queries, string name_doc_expanded_queries){

	
	vector<char> st1(max_size);
	vector< vector<char> > bestw;
	bestw.resize( max_size , vector<char>( N , 0 ) );
	vector< vector<char> > st;
	st.resize( max_size , vector<char>( 100 , 0 ) );
	float dist, len;
	vector<float> bestd(N);
	vector<float> vec(max_size);
	long long a , b, c, d, cn;
	long long words = vocab.size();
	cout<<"Size of vocab : " << words <<endl;
	vector<long long> bi(100);
	string line;
	int compteur = 0;
	

	ofstream doc_expanded_queries(name_doc_expanded_queries);
	ifstream doc_queries(name_doc_queries);
	
	while ( getline (doc_queries,line) ){
		compteur++;
    		for (long a = 0; a < N; a++) bestd[a] = 0;
    		for (long a = 0; a < N; a++) bestw[a][0] = 0;
		for(a = 0; a < int(line.size()) ; a++) st1[a] = line[a];
		st1[line.size()] = 0;
    		if (!strcmp(&st1[0], "EXIT")) break;
    		cn = 0;
    		b = 0;
    		c = 0;

	//Store each word in st1 in st
    		while (1) {
      			st[cn][b] = st1[c];
      			b++;
      			c++;
      			st[cn][b] = 0;
      			if (st1[c] == 0) break;
      			if (st1[c] == ' ') {

        			cn++;
       				b = 0;
        			c++;

      			}
    		}

    		cn++;

	//For each word in st check if it is in the vocabulary if a word is not in the vocabulary, the query will not be expanded
    		for (int a = 0; a < cn; a++) {

      			for (b = 0; b < words; b++) if (vocab[b]==&st[a][0]) break;
      			if (b == words) b = -1;
      			bi[a] = b;
      			printf("\nWord: %s Position in vocabulary: %lld\n", &st[a][0], bi[a]);
      			if (b == -1) {

        			printf("Out of dictionary word!\n");
				doc_expanded_queries << line + "\n";
        			break;

      			}
    		}

    		if (b == -1) continue;
    		printf("\n                                              Word       Cosine distance\n----------------------------------------------------------------\n");

	//If st is composed of 1 word, fill vec with the vectorial representation of st, otherwise, the sentence will be represented by the sum of the vector composing it 
    		for (uint64_t a = 0; a < size; a++) vec[a] = 0;
    		for (b = 0; b < cn; b++) {

      			if (bi[b] == -1) continue;
      			for (uint64_t a = 0; a < size; a++) vec[a] += M[a + bi[b] * size];

    		}

	//Normalization
    		len = 0;
    		for (uint64_t a = 0; a < size; a++) len += vec[a] * vec[a];
    		len = sqrt(len);
    		for (uint64_t a = 0; a < size; a++) vec[a] /= len;

    		for (a = 0; a < N; a++) bestd[a] = 0;
    		for (a = 0; a < N; a++) bestw[a][0] = 0;

	//For each word in the vocabulary 
    		for (c = 0; c < words; c++) {
			
		//Skip the identical terms 
      			a = 0;
	      		for (b = 0; b < cn; b++) if (bi[b] == c) a = 1;
	      		if (a == 1) continue;
	      		dist = 0;

		//Compute the distance with st and check if and where it is in the top N
	      		for (uint64_t a = 0; a < size; a++) dist += vec[a] * M[a + c * size];
	      		for (a = 0; a < N; a++) {

				if (dist > bestd[a]) {

		  			for (d = N - 1; d > a; d--) {

		    				bestd[d] = bestd[d - 1];
		    				strcpy(&bestw[d][0], &bestw[d - 1][0]);

		  			}

		  			bestd[a] = dist;
		  			strcpy(&bestw[a][0], &(vocab[c][0]));
		  			break;
				}
	      		}
    		}
	
	//Write in the file expanded_queries the new queries
    		for (a = 0; a < N; a++){ 

			printf("%50s\t\t%f\n", &bestw[a][0], bestd[a]);
			doc_expanded_queries << line + " " + string(&bestw[a][0]) + "\n";

  		}
  	}

	
	doc_queries.close();
	doc_expanded_queries.close();


}




void get_size(std::string file_name,size_t& words, size_t& size)
{
	
	FILE *f = fopen(file_name.c_str(), "rb");

	if(f == NULL) throw int(err1);

	
	fscanf(f, "%ld", &words);
  	fscanf(f, "%ld", &size);

	fclose(f);

}



template<typename T>
void ecriture_fichier_bin(T *var, string nom_fichier, long long size)
{


	FILE*f=fopen(nom_fichier.c_str(),"wb");

	fwrite(var ,sizeof(T),size,f);
	
	fclose(f);

}

template<typename T>
void lecture_fichier_bin(T *var, string nom_fichier, long long size)
{


	FILE *f = fopen(nom_fichier.c_str(), "rb");

	fread(var, sizeof(T), size, f);
	
	fclose(f);
	
}


void test_sentence(){

	FILE*f=fopen("true_queries.txt","rb");
	Sentence sentence;
	sentence.read_file(f);
	sentence.split();
	for(unsigned int i = 0; i < sentence.size() ; i++){

		cout<< i << " "<< sentence[i] <<endl;

	}

}

/*

int main(int argc, char **argv) {


	//test_sentence();
	//return 0;

	
//Warning message if there is no file 

  	if (argc < 2) {
    		printf("Usage: ./distance <FILE>\nwhere FILE contains word projections in the BINARY FORMAT\n");
    		return 0;
  	}

//Declaration
	// Name of the binary file containing the embeded vectors
	std::string file_name;
	file_name = std::string(argv[1]);

	Embedding embedding;
	cout<<"loading "<< file_name <<endl;

	clock_t begin = clock();

	embedding.load_Word2VecBinFormat(file_name);

	clock_t end = clock();
  	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<< "Time to do all the stuff : "<< elapsed_secs <<endl;
	
	embedding.display_attributes();

	cout<<"Cosine similarity between test and tests : "<<embedding.cosine("test","tests") <<endl;

	return 0;

	cout<<"Query : "<< argv[2] <<endl;
	cout<<"Exp query "<< argv[3] <<endl;
	embedding.expand_queries(argv[2], argv[3]);





//Reads the file in input and get the number of words and their associated vectors's length in the variables words and size

	// Vocabulary size
	size_t words = 0;
	// Embeded vectors size
	size_t size = 0;
	// Read in tghe header the two sizes
  	get_size( file_name, words , size);


//Allocation of memory for the variables vocab that will contain the vocabulary and M the matrix that will contain their associated vector

    // A basic word
    typedef char Word[max_w];
	// the table of all words
	vector<Word> vocab(words);

	// Matrix of all embeded vectors
	vector<float> M(words * size);
	
//Check if the memory of the matrix has been allocated 

  	if (M.size() == 0) {
    		printf("Cannot allocate memory: %lu MB    %lu  %lu\n", words * size * sizeof(float) / 1048576, words, size);
    		return -1;
  	}
	
//Read the vocabulary and the vectors 

	lecture_fichier_bin(&M[0], "vectors.bin" ,  size*words);

	lecture_fichier_bin(&vocab[0], "vocabulaire.bin" ,  max_w*words);

//Take a file of queries and write it's "expanded" version	

	//expand_queries("true_queries.txt", "expanded_queries.txt", &vocab,  &M, words, size);

  	return 0;
}

*/
