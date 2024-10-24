
#define DISCRIPTION_LENGTH     15 // constant definition (via preprocessor directive, old C feature)
#define NUM_NEURONS     2
unsigned long int myTime;
unsigned int mydelay = 10; // ms
/******************************************************/ 
/* Exercise 3, based on Mutually inhibiting neurons in Matsuoka 1985
 */
/******************************************************/ 
struct Neuron {
  double x = 0.0;
  double adaptation = 0.0;
  double adaptation_tau = 1.0;
  double y = 0.0;
  double b = 0.0; // scales the adaptation term
  double inj_cur = 0.0;

} neurons[NUM_NEURONS]; // we have control over individual neurons in the network


struct {
  double y[NUM_NEURONS];
  int connectionMatrix[NUM_NEURONS][NUM_NEURONS]  = {{0, 1}, {1, 0}};
} all_neurons;
/******************************************************/ 
// Non-linear function for y

double sigmoid(double x_i)
{
  return double (1/(1 + exp(-1*x_i)));
}

/******************************************************/ 
// dx/dt function
// synaptic strength must be submitted as a row vector that corresponds to coefficients
double fun_dx (double x_i, const int synaptic_strength[], const double y[], int len_y, double s_i , double b, double adaptation_i)
{ 
  double sum_y_scaled = 0.0;
  for(int j = 0; j < len_y; j++)
    {
    sum_y_scaled += y[j] * synaptic_strength[j];
    }
  return  (double)(-x_i - sum_y_scaled + s_i - b*adaptation_i); 
  } 


// d adaptation/dt function
double fun_adaptation (double a_i, double adaptation_tau, double y_i)
{ 
  return  (double)(1/adaptation_tau * (y_i - a_i)); 
  }   

/******************************************************/ 
/* put your setup code in setup(), to run once */
void setup() {
  Serial.begin(115200);
  for (int i =0; i < NUM_NEURONS; i++){
    all_neurons.y[i] = neurons[i].y;
  }

}
/******************************************************/
void update_one_neuron(struct Neuron* neuron_pointer, int neuron_index)
{
  int NUM_UPDATES = 2;
  double step; 
  step = ((double)mydelay/1000)/NUM_UPDATES;
  double x_i[NUM_UPDATES], adaptation_i[NUM_UPDATES]; 
  double k1,k2,k3,k4,k, l1,l2,l3,l4,l;  

  x_i[0] = neuron_pointer -> x; 
  adaptation_i[0] = neuron_pointer -> adaptation;

  int synaptic_s[NUM_NEURONS]; // only taking the row corresponding to the neuron we are updating
  for (int j = 0; j < NUM_NEURONS; j++) {
    synaptic_s[j] = all_neurons.connectionMatrix[neuron_index][j];
}



  // RUNGE-KUTTA ODE SOLVER
  for (int i=0; i< NUM_UPDATES-1 ; i++)
  {
    k1 = step*fun_dx(x_i[i], synaptic_s, all_neurons.y, NUM_NEURONS, neuron_pointer -> inj_cur, neuron_pointer -> b, adaptation_i[i]); //(double x_i, double synaptic_strength[], double y[], int len_y, double s_i , double b, double adaptation_i)
    l1 = step*fun_adaptation(adaptation_i[i], neuron_pointer -> adaptation_tau, neuron_pointer -> y); //(double adaptation_i, double adaptation_tau, y_i)

    k2 = step*fun_dx(x_i[i] + k1/2, synaptic_s, all_neurons.y,NUM_NEURONS, neuron_pointer -> inj_cur, neuron_pointer -> b, adaptation_i[i]+l1/2); //(double x_i, double synaptic_strength[], double y[], int len_y, double s_i , double b, double adaptation_i)
    l2 = step*fun_adaptation(adaptation_i[i] + l1/2, neuron_pointer -> adaptation_tau, neuron_pointer -> y); //(double adaptation_i, double adaptation_tau, y_i)

    k3 = step*fun_dx(x_i[i] + k2/2, synaptic_s, all_neurons.y, NUM_NEURONS, neuron_pointer -> inj_cur, neuron_pointer -> b, adaptation_i[i]+l2/2); //(double x_i, double synaptic_strength[], double y[], int len_y, double s_i , double b, double adaptation_i)
    l3 = step*fun_adaptation(adaptation_i[i]+l2/2, neuron_pointer -> adaptation_tau, neuron_pointer -> y); //(double adaptation_i, double adaptation_tau, y_i)

    k4 = step*fun_dx(x_i[i]+k3, synaptic_s, all_neurons.y, NUM_NEURONS, neuron_pointer -> inj_cur, neuron_pointer -> b, adaptation_i[i]+l3); //(double x_i, double synaptic_strength[], double y[], int len_y, double s_i , double b, double adaptation_i)
    l4 = step*fun_adaptation(adaptation_i[i]+l3, neuron_pointer -> adaptation_tau, neuron_pointer -> y); //(double adaptation_i, double adaptation_tau, y_i)   
    
    k = 1/6.0 * (k1 + 2*k2 + 2*k3 + k4);
    l  = 1/6.0 * ( l1 + 2*l2  + 2*l3  +  l4);
  
    x_i[i+1]= x_i[i]+k;
    adaptation_i[i+1] = adaptation_i[i]+l;
  } 

  neuron_pointer-> x = x_i[NUM_UPDATES-1]; 
  neuron_pointer-> adaptation = adaptation_i[NUM_UPDATES-1]; 
  double new_y  = sigmoid(x_i[NUM_UPDATES-1]);
  neuron_pointer-> y = new_y;
  all_neurons.y[neuron_index] = new_y;

  return; 

}

/******************************************************/ 
void update_network(void)
{
for (int i = 0; i< NUM_NEURONS ; i++)
  update_one_neuron(&neurons[i], i);
}

/******************************************************/ 
/* put your main code here in loop(), to run repeatedly */
void loop() {

  /* Read my program running time in milliseconds */
  myTime = millis();

  /* Update the neurons output*/
  for (int i = 0; i< NUM_NEURONS; i++) // give all  neurons current simulatneously
    {
      if((myTime>5000)&&(myTime<5200))
        {neurons[i].inj_cur = 1;}
      else
        {neurons[i].inj_cur = 0;}
    }


  update_network();


  /* Printing the output of the neurons on serial port*/
  for (int i = 0; i < NUM_NEURONS ; i++)
  {Serial.print(neurons[i].y);Serial.print(" ");}
  Serial.print("\n");

  /* delay at the end */
  delay(mydelay);

}
