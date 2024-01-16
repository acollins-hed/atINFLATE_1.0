void Initialize_Environment(){
  for(int i=0;i<S;i++){
    for(int j=0;j<S;j++){
      W(i,j) = pow(phi,abs(sitetypes(i)-sitetypes(j)));
    }
  }
  
  if(codon_ring_space){
    for(int i = 0;i<T;i++){
      for(int j = 0; j<T;j++){
	if(i == j)
	  mutation(i,j) = 1 - 2*Mu;
	else{
	  if(j == i+1 || j == i-1 || (j == 0 && i == T - 1) || (j == T - 1 && i == 0))
	    mutation(i,j) = Mu;
	  else
	    mutation(i,j) = 0;
	}
      }
    }
  }
  else{
    if(nbase == 2){
      Eigen::MatrixXd Codon_Mutation(2,2);
      for(int i=0;i<2;i++){
	for(int j=0;j<2;j++){
	  if(i == j)
	    Codon_Mutation(i,j) = 1-Mu;
	  else
	    Codon_Mutation(i,j) = Mu;
	}
      }

      if(T == 2)
	mutation = Codon_Mutation;
      if(T == 4){
	mutation = kroneckerProduct(Codon_Mutation,Codon_Mutation);
      }
      if(T == 8){
	mutation = kroneckerProduct(Codon_Mutation,kroneckerProduct(Codon_Mutation,Codon_Mutation));
      }
    }
    if(nbase == 4){
      Eigen::MatrixXd Codon_Mutation(4,4);
      
      for(int i=0;i<4;i++){
	for(int j=0;j<4;j++){
	  if(i == j)
	    Codon_Mutation(i,j) = 1-Mu;
	  else{
	    if((i == 0 && j == 1)||(i == 1 && j==0) ||(i == 2 && j == 3)||(i == 3 && j == 2))
	      Codon_Mutation(i,j) = ((double) kappa) * Mu/((double)kappa+2);
	    else
	      Codon_Mutation(i,j) = Mu/((double) kappa + 2);
	  }
	}
      }
      
      if(T == 4)
	mutation = Codon_Mutation;
      else{
	
	if(T == 16){
	  
	  mutation = kroneckerProduct(Codon_Mutation,Codon_Mutation);
	  
	}
	if(T == 64)
	  mutation = kroneckerProduct(Codon_Mutation,kroneckerProduct(Codon_Mutation,Codon_Mutation));
      }      
    }
  }
  cout<<"The Codon Mutations are\n"<<mutation<<endl;
}
		     
double hmean(Eigen::MatrixXd k, int i){
  double sum=0;
  for(int j=0;j<Acap;j++){
	sum += 1/k.coeff(i,j);
  }
  return (double(Tcap)/(sum));
}

double A_fit(Eigen::MatrixXd* codon_freq, Eigen::MatrixXd* eff_mutant_code, int stype){
  double f=0;

  for(int aa=0;aa<Acap;aa++){
    for(int codon=0;codon<Tcap;codon++){
      f += W(stype,aa_to_st[amino_acids(aa)])*eff_mutant_code->coeff(codon,aa)*codon_freq->coeff(stype,codon);
    }
  }
  return f;
}

double M_fit(Eigen::MatrixXd* codon_freq, Eigen::MatrixXd* eff_mutant_code, Eigen::MatrixXd* kd){
  double f=1;

  for(int stype=0;stype<S;stype++){
    f *= pow(A_fit(codon_freq,eff_mutant_code,stype),l(stype));
  }

  if(rate){
    double rexponent = 1.505764/44;
    double kdbar=0;

    for(int i=0; i<Tcap;i++){
      for(int j=0; j<Acap; j++){
	kdbar+=1/kd->coeff(i,j);
      }
    }
    kdbar = Tcap*Acap/kdbar;
    
    //////////////////
    
    f *= (1 - exp(-rexponent*kdbar));
  }
  return f;
}

struct tRNAs
{
  Eigen::VectorXi anticodons;
  //The iis of the tRNAs are 2xT matrices.
  Eigen::MatrixXi iis;
  tRNAs(){
    Eigen::MatrixXi t(2,T);
    Eigen::VectorXi a(T);
    t.setZero();
    iis = t;
    for(int i=0;i<T;i++)
      a(i) = i;
    anticodons=a;
  }
  tRNAs(Eigen::MatrixXi x, Eigen::VectorXi z){
    iis = x;
    anticodons = z;
  }
  void  print(){
    for(int i = 0;i<2;i++){
      for(int j=0;j<T;j++){
	for(int ii=n-1;ii>=0;ii--){
	  if(iis(i,j)&(1<<ii))
	    cout<<"1";
	  else
	    cout<<"0";
	}
	cout<<" ";
      }
      cout<<endl;
    }
  }

  string print_trna(int i, int j){
    string s = "";
    for(int ii=n-1;ii>=0;ii--){
      if(iis(i,j)&(1<<ii))
	s += "1";
      else
	s += "0";
    }
    return s;
  }
  
};

struct aaRSs
{
  Eigen::VectorXf aas;
  //The iis of the aaRSs are 2xA matrices.
  Eigen::MatrixXi iis;
  aaRSs(){
    Eigen::MatrixXi a(2,A);
    Eigen::VectorXf aa(A);
    a.setZero();
    iis = a;
    for(int i=0;i<A;i++)
      aa(i) = ((float)i)/(A-1);
    aas=aa;
  }
  aaRSs(Eigen::MatrixXi x, Eigen::VectorXf y){
    iis = x;
    aas = y;
  }
  void  print(){
    for(int i = 0;i<2;i++){
      for(int j=0;j<A;j++){
	for(int ii=n-1;ii>=0;ii--){
	  if(iis(i,j)&(1<<ii))
	    cout<<"1";
	  else
	    cout<<"0";
	}
	cout<<" ";
      }
      cout<<endl;
    }
  }

  string print_aars(int i, int j){
    string s = "";
    for(int ii=n-1;ii>=0;ii--){
      if(iis(i,j)&(1<<ii))
	s += "1";
      else
	s += "0";
    }
    return s;
  }

};

struct Current_Genotype
{
 public:
  tRNAs trnas;
  aaRSs aarss;
  Eigen::MatrixXd code;
  Eigen::MatrixXd kd;
  void Get_Code(){
    kd=Get_kd();
    Eigen::MatrixXd c(T,A);
    double hx;
    for(int i=0;i<Tcap;i++){
      hx = hmean(kd,i)/double(Tcap);
      for(int j=0;j<Acap;j++){
	c(i,j) = hx/kd(i,j);
      }
    }
    code = c;
  }
  
  Current_Genotype(){
    tRNAs t;
    aaRSs a;
    trnas=t;
    aarss=a;
    Get_Code();
  }
  Current_Genotype(tRNAs t, aaRSs a)
  :trnas(t),aarss(a){
    Get_Code();
  }
  Eigen::MatrixXd Get_kd(){
    Eigen::MatrixXd krate(T,A);
    krate.setZero();

    for(int i=0;i<Tcap;i++){
      for(int j=0;j<Acap;j++){
	if(__builtin_popcount((trnas.iis(1,i)&aarss.iis(1,j))&(~(trnas.iis(0,i)^aarss.iis(0,j)))) > k)
	  krate(i,j) = kdc;
	else
	  krate(i,j) = kdnc*exp(-1*epsilon*(__builtin_popcount((trnas.iis(1,i)&aarss.iis(1,j))&(~(trnas.iis(0,i)^aarss.iis(0,j))))));
      }
    }

    return krate;
  }
};

struct Population
{
  Current_Genotype genotype;
  double fitness;
  Eigen::MatrixXd codon_frequency;
  Population(){
    Current_Genotype genotype;
    Get_Codon_Freq();
    fitness = M_fit(&codon_frequency, &genotype.code, &genotype.kd);
  }
  Population(Current_Genotype g)
    :genotype(g){
    Get_Codon_Freq();
    fitness = M_fit(&codon_frequency, &genotype.code, &genotype.kd);
  }
  void Get_Codon_Freq(){
    Eigen::MatrixXd Q(T,T),wfit(T,T),diag(S,T);
    diag = Diagonal();
    codon_frequency.resize(S,T);
    for(int i=0;i<S;i++){
      wfit.setZero();
      wfit.diagonal() = diag.row(i).array();
      Q = mutation*wfit;
      Eigen::EigenSolver<Eigen::MatrixXd> ev(Q);
      //I have seen claims on the internet that Eigen automatically puts the
      //largest eigenvalue into cell 0 (I've not seen this on Eigen's own
      //website) and I've found this to be untrue, hence the following line.
      codon_frequency.row(i) = ev.eigenvectors().col(max_eigenvalue(&ev.eigenvalues())).real();
      //The following makes sure the components sum to one. The eigenvectors are
      //already normalized but of course that just makes their norm 1, not their
      //components sum to one.
      Sum_to_one(&codon_frequency);
    }
  }

private:
  void Sum_to_one(Eigen::MatrixXd * emat){
    double sum=0;
    for(int i=0;i<emat->rows();i++){
      for(int j=0;j<emat->cols();j++){
	sum = emat->coeff(i,j) + sum;
      }
      emat->row(i) *= (1/sum);
	sum = 0;
    }
  }

  Eigen::MatrixXd Diagonal(){
    Eigen::MatrixXd diag(S,T);
    for(int stype = 0;stype<S;stype++){
      for(int codon=0;codon<Tcap;codon++){
	diag(stype,codon) = 0;
	for(int alpha=0;alpha<Acap;alpha++){
	  diag(stype,codon) += W(stype,aa_to_st[amino_acids(alpha)])*genotype.code(codon,alpha);
	}
      }
    }
    return diag;
  }

  int max_eigenvalue(const Eigen::EigenSolver<Eigen::MatrixXcd>::EigenvalueType * evals){
    int max = 0,i;
    for(i=1;i<evals->size();i++){
      if(evals->coeff(max).real() < evals->coeff(i).real()){
	max = i;
      }
    }
    return max;
  }
  
};

double Transition_Probability(double current_fitness, double mutant_fitness,int hamming){
  double Fixation_Probability;
  if(current_fitness == mutant_fitness)
    Fixation_Probability = 1/((double) N);
  else
    Fixation_Probability = (1 - (current_fitness)/(mutant_fitness))/(1 - pow((current_fitness)/(mutant_fitness),N));
  return (N*pow(mu,hamming)*Fixation_Probability); 
}

struct Evolver
{
  Eigen::MatrixXd single_mutant, double_mutant_wout, double_mutant_win;
  Population population;
  double current_fitness;

  Evolver(Population p)
    :population(p){
    current_fitness = population.fitness;
    single_mutant = Eigen::MatrixXd(L,5);
    double_mutant_win = Eigen::MatrixXd(n_d_win,5);
    double_mutant_wout = Eigen::MatrixXd(n_d_wout,3);
    Get_Mutant_Web();
  }

  void Get_Mutant_Web(){
    double sum = 0;
    int index = 0;
    single_mutant.setZero();
    double_mutant_wout.setZero();
    double_mutant_win.setZero();
    
    //Beginning of single_mutant loops
    //
    //column 0 = 0 for tRNA, 1 for aaRS
    //column 1 = row of mutation
    //column 2 = column of mutation
    //column 3 = value in (row,column)
    //column 4 = transition probability

    for(int row=0;row<2;row++){
      for(int col=0;col<Tcap;col++){
	for(int i=0;i<n;i++){
	  Population mutant_pop = population;
	  mutant_pop.genotype.trnas.iis(row,col) ^= 1<<i;
	  mutant_pop.genotype.Get_Code();
	  mutant_pop.fitness = M_fit(&population.codon_frequency,&mutant_pop.genotype.code,&mutant_pop.genotype.kd);
	  single_mutant(index,0) = 0;
	  single_mutant(index,1) = row;
	  single_mutant(index,2) = col;
	  single_mutant(index,3) = mutant_pop.genotype.trnas.iis(row,col);
	  single_mutant(index,4) = Transition_Probability(current_fitness,mutant_pop.fitness,1);
	  sum += single_mutant(index,4);
	  index++;
	}
      }
    }

    for(int row=0;row<2;row++){
      for(int col=0;col<Acap;col++){
	for(int i=0;i<n;i++){
	  Population mutant_pop = population;
	  mutant_pop.genotype.aarss.iis(row,col) ^= 1<<i;
	  mutant_pop.genotype.Get_Code();
	  mutant_pop.fitness = M_fit(&population.codon_frequency,&mutant_pop.genotype.code,&mutant_pop.genotype.kd);
	  single_mutant(index,0) = 1;
	  single_mutant(index,1) = row;
	  single_mutant(index,2) = col;
	  single_mutant(index,3) = mutant_pop.genotype.aarss.iis(row,col);
	  single_mutant(index,4) = Transition_Probability(current_fitness,mutant_pop.fitness,1);
	  sum += single_mutant(index,4);
	  index++;
	}
      }
    }
    
    //Beginning of double_mutant_win loops
    //
    //column 0 = 0 for tRNA, 1 for aaRS
    //column 1 = row of mutation
    //column 2 = column of mutation
    //column 3 = value in (row,column)
    //column 4 = transition probability
    index = 0;
    for(int row=0;row<2;row++){
      for(int col=0;col<Tcap;col++){
	for(int i=0;i<n-1;i++){
	  for(int j=n-1;j>i;j--){
	    Population mutant_pop = population;
	    mutant_pop.genotype.trnas.iis(row,col) ^= 1<<i;
	    mutant_pop.genotype.trnas.iis(row,col) ^= 1<<j;
	    mutant_pop.genotype.Get_Code();
	    mutant_pop.fitness = M_fit(&population.codon_frequency,&mutant_pop.genotype.code,&mutant_pop.genotype.kd);
	    double_mutant_win(index,0) = 0;
	    double_mutant_win(index,1) = row;
	    double_mutant_win(index,2) = col;
	    double_mutant_win(index,3) = mutant_pop.genotype.trnas.iis(row,col);
	    double_mutant_win(index,4) = Transition_Probability(current_fitness,mutant_pop.fitness,2);
	    sum += double_mutant_win(index,4);
	    index++;
	  }
	}
      }
    }

    for(int row=0;row<2;row++){
      for(int col=0;col<Acap;col++){
	for(int i=0;i<n-1;i++){
	  for(int j=n-1;j>i;j--){
	    Population mutant_pop = population;
	    mutant_pop.genotype.aarss.iis(row,col) ^= 1<<i;
	    mutant_pop.genotype.aarss.iis(row,col) ^= 1<<j;
	    mutant_pop.genotype.Get_Code();
	    mutant_pop.fitness = M_fit(&population.codon_frequency,&mutant_pop.genotype.code,&mutant_pop.genotype.kd);
	    double_mutant_win(index,0) = 1;
	    double_mutant_win(index,1) = row;
	    double_mutant_win(index,2) = col;
	    double_mutant_win(index,3) = mutant_pop.genotype.aarss.iis(row,col);
	    double_mutant_win(index,4) = Transition_Probability(current_fitness,mutant_pop.fitness,2);
	    sum += double_mutant_win(index,4);
	    index++;
	  }
	}
      }
    }


    //Beginning of double_mutant_wout loops
    //
    //column 0 = row of first mutation in single_mutant
    //column 1 = row of second mutation in single_mutant
    //column 2 = transition probability

    //An explanation for the indexing of the following set of nested
    //loops. i and j are used to pull mutation information from
    //single_mutant. i begins at 0 and ticks up to L-n-1. j begins
    //at q*n+n and ticks up to L-1 where  q is the quotient i with n, i.e.
    //i = qn + r. This way, j begins at n bit multiples and is at least
    //one interaction interface ahead of i.
    index = 0;

    for(int i=0;i<L-n;i++){
      for(int j = n*(i/n)+n;j<L;j++){
	double_mutant_wout(index,0) = i;
	double_mutant_wout(index,1) = j;
	Population mutant_pop = population;
	if(single_mutant(i,0))
	  mutant_pop.genotype.aarss.iis(single_mutant(i,1),single_mutant(i,2)) = single_mutant(i,3);
	else
	  mutant_pop.genotype.trnas.iis(single_mutant(i,1),single_mutant(i,2)) = single_mutant(i,3);
	if(single_mutant(j,0))
	  mutant_pop.genotype.aarss.iis(single_mutant(j,1),single_mutant(j,2)) = single_mutant(j,3);
	else
	  mutant_pop.genotype.trnas.iis(single_mutant(j,1),single_mutant(j,2)) = single_mutant(j,3);
	mutant_pop.genotype.Get_Code();
	mutant_pop.fitness = M_fit(&population.codon_frequency,&mutant_pop.genotype.code,&mutant_pop.genotype.kd);
	double_mutant_wout(index,2) = Transition_Probability(current_fitness,mutant_pop.fitness,2);
	sum += double_mutant_wout(index,2);
	index++;
      }
    }

    single_mutant.col(4) = (1/sum)*single_mutant.col(4);
    double_mutant_win.col(4) = (1/sum)*double_mutant_win.col(4);
    double_mutant_wout.col(2) = (1/sum)*double_mutant_wout.col(2);
  }

  void Evolve(){
    double rand_n = dist(mt), sum=0;
    bool transition = 0;

    Get_Mutant_Web();
    //Checking if the single mutants transition probabilities
    //exceed the random number
    
    for(int i=0;i<L;i++){
      sum+=single_mutant(i,4);
      if(sum > rand_n){
	if(single_mutant(i,0))
	  population.genotype.aarss.iis(single_mutant(i,1),single_mutant(i,2)) = single_mutant(i,3);
	else
	  population.genotype.trnas.iis(single_mutant(i,1),single_mutant(i,2)) = single_mutant(i,3);
	population.genotype.Get_Code();
	population.Get_Codon_Freq();
	population.fitness = M_fit(&population.codon_frequency,&population.genotype.code,&population.genotype.kd);
	transition = 1;
	mutation_type="Single";
	break;
      }
    }

    //Checking if the double mutants_win transition probabilities
    //exceed the random number
    if(!transition){
      for(int i=0;i<n_d_win;i++){
	sum+=double_mutant_win(i,4);
	if(sum > rand_n){
	  if(double_mutant_win(i,0))
	    population.genotype.aarss.iis(double_mutant_win(i,1),double_mutant_win(i,2)) = double_mutant_win(i,3);
	  else
	    population.genotype.trnas.iis(double_mutant_win(i,1),double_mutant_win(i,2)) = double_mutant_win(i,3);
	  population.genotype.Get_Code();
	  population.Get_Codon_Freq();
	  population.fitness = M_fit(&population.codon_frequency,&population.genotype.code,&population.genotype.kd);
	  transition = 1;
	  mutation_type="Double_win";
	  break;
	}
      }
    }

    //Checking if the double mutants_wout transition probabilities
    //exceed the random number
    if(!transition){
      for(int i = 0;i<n_d_wout;i++){
	sum+=double_mutant_wout(i,2);
	if(sum>rand_n){
	  if(single_mutant(double_mutant_wout(i,0),0))
	    population.genotype.aarss.iis(single_mutant(double_mutant_wout(i,0),1),single_mutant(double_mutant_wout(i,0),2)) = single_mutant(double_mutant_wout(i,0),3);
	  else
	    population.genotype.trnas.iis(single_mutant(double_mutant_wout(i,0),1),single_mutant(double_mutant_wout(i,0),2)) = single_mutant(double_mutant_wout(i,0),3);
	  if(single_mutant(double_mutant_wout(i,1),0))
	    population.genotype.aarss.iis(single_mutant(double_mutant_wout(i,1),1),single_mutant(double_mutant_wout(i,1),2)) = single_mutant(double_mutant_wout(i,1),3);
	  else
	    population.genotype.trnas.iis(single_mutant(double_mutant_wout(i,1),1),single_mutant(double_mutant_wout(i,1),2)) = single_mutant(double_mutant_wout(i,1),3);
	  population.genotype.Get_Code();
	  population.Get_Codon_Freq();
	  population.fitness = M_fit(&population.codon_frequency,&population.genotype.code,&population.genotype.kd);
	  transition = 1;
	  mutation_type="Double_out";
	  break;
	}
      }
    }
    population.genotype.Get_Code();
    population.Get_Codon_Freq();
    population.fitness = M_fit(&population.codon_frequency,&population.genotype.code,&population.genotype.kd);
    current_fitness = population.fitness;
  }
  
};

