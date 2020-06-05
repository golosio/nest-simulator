import nest
import math
import matplotlib.pyplot as plt

sim_time = 20.0

# STDP weight update theoretical formula for comparison
def STDPUpdate(w, Dt, tau_plus, tau_minus, Wplus, alpha, mu_plus, mu_minus, \
               Wmax):
    if (Dt>=0):
        fact = Wplus*math.exp(-Dt/tau_plus)
        w1 = w + fact*math.pow(1.0 - w/Wmax, mu_plus)
        if w1>Wmax:
            w1 = Wmax
        
    else:
        fact = -alpha*Wplus*math.exp(Dt/tau_minus)
        w1 = w + fact*math.pow(w/Wmax, mu_minus)
        if w1<0.0:
            w1 = 0.0
    return w1


# presynaptic and postsynaptic neurons
neuron_pre = nest.Create("parrot_neuron")
neuron_post = nest.Create("parrot_neuron")

#spike generators
sg_pre = nest.Create("spike_generator")
sg_post = nest.Create("spike_generator")


# spike times
spike_times_pre = [4.0]
spike_times_post = [1.0, 5.0]
nest.SetStatus(sg_pre, {"spike_times": spike_times_pre})
nest.SetStatus(sg_post, {"spike_times": spike_times_post})

# connect spike generators to neurons
syn_dict={"weight":1.0, "delay":1.0}
nest.Connect(sg_pre, neuron_pre, "one_to_one", syn_dict)
nest.Connect(sg_post, neuron_post, "one_to_one", syn_dict)

# STDP connection parameters
tau_plus = 20.0
tau_minus = 20.0
lambd = 0.01
alpha = 1.0
mu_plus = 1.0
mu_minus = 1.0
Wmax = 10.0
den_delay = 3.0 
weight_stdp = 1.0

delay = 1.0
syn_dict_stdp={"weight":weight_stdp, "delay":delay, \
               "synapse_model":"stdp_synapse",
               "tau_plus":tau_plus, \
               "lambda":lambd, "alpha":alpha, "mu_plus":mu_plus, \
               "mu_minus":mu_minus,  "Wmax":Wmax, "den_delay":den_delay,
               # set receptor 1 postsynaptically, to not generate extra spikes
               "receptor_type": 1}

nest.Connect(neuron_pre, neuron_post, "one_to_one", syn_dict_stdp)

nest.Simulate(sim_time)

conn_id = nest.GetConnections(neuron_pre, neuron_post)
w = nest.GetStatus(conn_id, "weight")

print("Initial weight: ", weight_stdp)
print("Simulated weight: ", w[0])

Wplus = Wmax*lambd
Dt1 = -1.0
w1 = STDPUpdate(weight_stdp, Dt1, tau_plus, tau_minus, Wplus, \
                alpha, mu_plus, mu_minus, Wmax)
Dt2 = 3.0
w2 = STDPUpdate(w1, Dt2, tau_plus, tau_minus, Wplus, \
                alpha, mu_plus, mu_minus, Wmax)

print("Expected theoretical weight: ", w2)

print("dw/w: ", (w2 - w[0])/w2)
