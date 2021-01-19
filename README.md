# Agent based simulation of the spread of COVID-19 in confined spaces
**Author: Jana Lasser, Complexity Science Hub Vienna (lasser@csh.ac.at)**

A simulation tool to explore the spread of COVID-19 in small communities such as nursing homes or schools via agent-based modeling (ABM) and the impact of prevention measures. The model follows an SEIRX approach, building on the agent based simulation framework [mesa](https://mesa.readthedocs.io/en/master/) in which agents can be susceptible (S), exposed (E), infected (I), removed (R) or quarantined (X) and is based on explicitly defined and dynamic contact networks between agents. The model offers the possibility to explore the effectiveness of various testing, tracing and quarantine strategies and other interventions such as ventilation and mask-wearing.

<img alt="Illustrative figure of infection spread in a nursing home" src="img/fig.png?raw=true" height="500" width="800" align="center">


**This software is under development and intended to respond rapidly to the current situation. Please use it with caution and bear in mind that there might be bugs**

Reference:  

_Lasser, J. (2020). Agent based simulation of the spread of COVID-19 in nursing homes. [DOI](https://doi.org/10.5281/zenodo.4275533): 10.5281/zenodo.4275533_

## Simulation design

### Agents
#### States
Simulations can have several types of agents, for example residents & employees in nursing homes, or teachers, students and family members in schools. Infections are introduced through agents that have a certain probability to become an index case or that are explicitly chosen as index case in the beginning of the simulation. The contact network defines which agents interact with which other agents on which day of the week. Different contact types modulate the infection transmission risk between "close" and "very far". Contact networks are stored as a [networkx](https://networkx.org/) graph with edge attributes for different contact types. In every step (day) of the simulation, agents interact according to the contact network and can transmit the infection to the agents they are in contact with. Depending on its infection state, an agent has one of four states: susceptible (S), exposed (E), infected (I) or removed (R). In addition, agents can be quarantined (X) and develop symptoms when they are infected. While infected (I), agents can by symptomatic or asymptomatic. 

#### Attributes
Next to the dynamic states pertaining to the infection state of an agent, agents also have (static) attributes that influence how they interact with the infection dynamics:
* Exposure duration: The time between transmission and becoming infectious (exposure duration) for every agent is drawn from a Weibull distribution and might be different for every agent (see section [Parameters](#parameters)) for details).
* Time until symptoms: The time between transmission and (potentially) showing symptoms for every agent is drawn from a Weibull distribution and might be different for every agent (see section [Parameters](#parameters) for details).
* Infection duration: The time between transmission and ceasing to be infectious is drawn from a Weibull distribution and might be different for every agent (see section [Parameters](#parameters) for details).
* Age: especially in the school setting, where children are involved, the age of the agents plays an important role, since children have a somewhat reduced risk to transmit or receive an infection (see section [Transmissions](#transmissions) for details).

#### Viral Load
An infection with SARS-CoV-2 causes an infected person's body to replicate the virus and the viral load in the person will vary, depending on the progression of the infection. In this simulation, viral load influences two distinct processes: the ability of testing technologies to detect an infection (detection threshold, see section [Test Technologies](#test-technologies)) and the infectiousness of agents (see section [Transmissions](#transmissions)). So far, we do not model the dependence of viral load on time explicitly for every agent. As pertaining to infectiousness, we model the viral load as trapezoid function that is high at the beginning and drops to zero over the course of the infection. For different testing technologies, we model the time between transmission and detection threshold individually for every testing technology used.

### Transmissions
Transmissions are modeled as [Bernoulli-Trials](https://en.wikipedia.org/wiki/Bernoulli_trial) with a probability p of success and a probability q = 1 - p of failure. In every step of the simulation (one step corresponds to one day), every infected and non-quarantined agent performs this Bernoulli trial once for every other (non-infected and non-quarantined) agent they are in contact with. The overall probability of a successful transmission between two agents is reduced by a range of factors that reflect both biological factors and intervention measures that can act on both the _transmitter_ of the infection as well as the _receiver_ of the infection. If q_m is the probability of failure of transmission due to wearing a mask and b is the baseline transmission risk without any modifications, then the modified probability of successful transmission between two agents is

p = 1 - [1 - b(1 - q_m)].

We currently account for eight different factors that can influence the transmission risk in different settings (see section [Parameters](#parameters) for details on how these factors can be set and how values for them are chosen):
* q_1: Modification of the transmission risk due to the type of contact between agents. Here, q_1 = 1 for a household contact (contact type "close", no modification), whereas q_1 < 1 for other contact types ("intermediate", "far", "very far", reduction of transmission risk). The value of q_1 depending on the type of the contact has to be specified via the model parameter ```infection_risk_contact_type_weights``` at model setup. Contact types between each two agents are stored in the contact network supplied to the model.
* q_2: Modification of the transmission risk due to the age of the _transmitting agent_. The dependence of transmission risk on age is set via the model parameter ```age_symptom_discount``` at model setup.
* q_3: Modification of the reception risk due to the age of the _receiving agent_. The dependence of reception risk on age is approximated to be the same as the dependence of the transmission risk on age and is therefore also set via the model parameter ```age_symptom_discount``` at model setup.
* q_4: Modification of the transmission risk due to the progression of the infection. This dependence is currently hard-coded, based on literature values.
* q_5: Modification of the transmission risk due to the type of the course of the infection (symptomatic, asymptomatic). This parameter is set via the model parameter ```subclinical_modifier``` at model setup.
* q_6: Modification of the transmission risk due to mask wearing of the _transmitting agent_. This parameter is set via the model parameter ```mask_filter_efficiency["exhale"]``` at model setup. Whether or not an agent group is wearing masks has to be specified via the ```agent_types["mask"]``` parameter and the contact types that are affected by mask-wearing are hard-coded in the model (for example, household contacts are generally not affected by mask-wearing).
* q_7: Modification of the reception risk due to mask wearing of the _receiving agent_. This parameter is set via the model parameter ```mask_filter_efficiency["inhale"]``` at model setup. Whether or not an agent group is wearing masks has to be specified via the ```agent_types["mask"]``` parameter and the contact types that are affected by mask-wearing are hard-coded in the model (for example, household contacts are generally not affected by mask-wearing).
* q_8: Modification of the transmission risk due to room ventilation. This parameter is set via the model parameter ```transmission_risk_ventilation_modifier``` at model setup.  

The baseline transmission risk is set via the model parameter ```base_transmission_risk``` at model setup.

Therefore, for example in a school setting where agents are wearing masks, rooms are ventilated and the age of agents is important for the transmission dynamics, the overall success probability for a transmission is defined as  

p = 1 - [1 - b(1 - q_1)(1 - q_2)(1-q_3)(1 - q_4)(1 - q_5)(1 - q_6)(1 - q_7)(1 - q_8)].

### Containment strategies
#### Testing strategies
Next to the transmission of the infection, containment measures (quarantine) and a testing and tracing strategy can be implemented to curb the spread of the virus. The general testing strategy can be specified via the model parameter ```testing```. 
* If ```testing='diagnostic'```, symptomatic cases are immediately quarantined and tested. Once a positive test result is returned, all K1 contacts of the positive agent are immediately quarantined too.  
* If ```testing='background'```, in addition to testing of single symptomatic agents, if there is a positive test result, a "background screen" of the population will be launched. In the nursing home scenario, all residents and employees are tested in such a background screen. In the school scenario, all teachers and students (but not family members) are tested in such a background screen. If a ```follow_up_testing```-interval is specified, each background screen is followed by a "follow up screen" that is similar to the background screen with the specified time-delay, testing the same agent groups as in the background screen and using the same testing technology.
* Next to population screening that is triggered by positive test results, if ```testing='preventive'```, preventive screens will be performed independently of diagnostic testing and background/follow-up screens. These preventive screens are performed in given intervals, which are to be specified for each agent group using the parameter ```agent_types[screening_interval]```. These screening intervals are tied to specific days of the week:
    * An interval of 7 days will cause preventive screens to be launched on Mondays
    * An interval of 3 days will cause preventive screens to be launched on Mondays and Thursdays
    * An interval of 2 days will cause preventive screens to be launched on Mondays, Wednesdays and Fridays.
    * An interval of None will not initiate any preventive screens, even if ```testing=preventive``` and will fall back to diagnostic testing and background screens.
    * Other intervals are currently not supported.
    
Tests take a certain amount of time to return results, depending on the chosen testing technology. Agents can have a pending test result, which will prevent them from getting tested again before the pending result arrives.  Tests can return positive or negative results, depending on whether the agent was testable at the time of testing (see section [Viral Load](#viral-load)) and on the sensitivity/specificity of the chosen test (see section [Test Technologies](#test-technologies)). 

#### Tracing
If an agent receives a positive test result (after the specified turnover time of a test, see [Test Technologies](#test-technologies)), their contacts are traced and also quarantined. The types of contacts between agents that will considered to be [K1](https://ehs.pages.ist.ac.at/definitions/) and will cause contact persons to be quarantined in case of a positive result of one of their contacts can be specified by the model parameter ```K1_contact_types```. Tracing is considered to occur instantly and contact persons are quarantined without time delay, as soon as a positive test result returns.

#### Quarantine
Quarantined agents will stay in quarantine for a number of days specified by the model parameter ```quarantine_duration``` days if ```liberating_testing=False``` (default). If ```liberating_testing = True```, quarantined agents will be released from quarantine if they return a negative test result. This has to be used with caution, as with test turnover times > 0 days, agents can have pending tests at the time they are quarantined, and a negative test result the next day or the day after can cause these agents to terminate their quarantine, even though they did not receive a test while in quarantine.

#### Test technologies
Depending on this progression, agents can be testable (i.e. tested positive) by different testing technologies, depending on the technology's detection threshold. In general, PCR tests are considered to be the gold-standard here, being able to detect very small viral loads (see section [Viral Load](#viral-load)), whereas antigen tests need considerably larger viral loads to detect an infection. In addition, tests can have a sensitivity and specificity, determining the probability to truthfully detect an infection (sensitivity) and the probability to truthfully determine a non-infection (specificity). Lastly, different test technologies need a different amount of time to return results. Here, antigen tests lead the field by only taking minutes to yield a result, whereas PCR tests require complex laboratory processing and can take several days until a result is found and communicated. 

A range of different test technologies such as ```same_day_antigen``` or ```two_day_PCR``` are specified in the file ```testing_strategy.py```. A test technology always specifies the test's sensitivity and specificity, the time until an agent is testable (from the day of transmission), the time an agent remains testable (from the day of transmission) and the test result turnover time.  

Test technologies for preventive screening and diagnostic testing (diagnostic tests, background screens and follow-up screens) can be specified separately, using the model parameters ```diagnostic_test_type``` and ```preventive_screening_test_type```.


## Implementation
### Model
The simulation consists of a _model_ that stores model parameters, the agent contact network and references to all agents. In every step of the simulation, the model initiates agent interactions and executes the testing and tracing strategy as well as data collection. Model parameters and parameters for the testing strategy, as well as a specification of the agent types and their respective parameters have to be passed to the  model instance at time of creation, if values other than the specified default values (see section [Parameters](#parameters)) should be used. The model is implemented as a class (```model_SEIRX```) that inherits from [mesa's](https://mesa.readthedocs.io/en/stable/) ```Model``` and implements a ```step()``` function that is called in every time-step. Every scenario (so far: nursing homes and schools) implements its own model class which inherits from ```model_SEIRX```, where functionality deviating from the behaviour specified in ```model_SEIRX``` is implemented. This can for example a custom ```step()``` or ```calculate_transmission_probability()``` function and specify scenario-specific data collection functions.

### Agents
Similarly to the model, agents have a base-class ```agent_SEIRX``` that inherits from [mesa's](https://mesa.readthedocs.io/en/stable/) ```Agent```. The agent class implements agent states and counters and functions necessary for simulating contacts between agents and advancing states. Different agent types needed in the scenarios inherit from this base class and might implement additional functionality. Currently there are five agent types: resident and employee (nursing home scenario) and teacher, student and family member (school scenario). These agent types are implemented in separate classes which inherit from the agent base-class. 

### testing
The testing strategy is contained in ```testing_strategy.py```, a class different from the SEIRX base model but is created with parameters passed through the SEIRX constructor. This is to keep parameters and information related to testing and tracing in one place, separate from the infection dynamics model. The testing class also stores information on the sensitivity, specificity and turnover time of a range of tests and can be easily extended to include additional testing technologies.

### Additional modules
* The module ```analysis_functions.py``` provides a range of functions to analyse data from model runs.
* The module ```viz.py``` provides some custom visualization utility to plot infection time-lines and agent states on a network, given a model instance.

## Applications
### Nursing homes
Nursing homes implement agents types ```resident``` and ```employee```, as well as the ```model_nursing_home``` (all located in the ```nursing_home``` sub-folder).  
The contact networks for nursing homes are specified through resident relations in the homes (room neighbors, table neighbors at joint meals, other shared areas). Employees interact with all other employees and all residents in the same living quarter of the nursing home every day. We provide several exemplary contact networks (see ```data/nursing_homes```), representing different architectures of nursing homes with different numbers of living quarters. These contact networks are abstractions of empirically determined interaction relations in Austrian nursing homes. By default, resident roommates are defined as "close contact", residents that share tables during joint meals as well as all resident - employee contacts are considered "intermediate", and contacts between residents that only share the same living quarters but not the same room or table are considered to be "far". 

### Schools
Schools implement agent types ```teachers```, ```students``` and ```family_members``` of students, as well as the ```model_school``` (all located in the ```school``` sub-folder).  

The contact networks for schools are modeled to reflect common structures in Austrian schools (see the [school-type specific documentation](https://github.com/JanaLasser/agent_based_COVID_SEIRX/blob/master/school/school_type_documentation.ipynb) for details). For this application, we construct three distinct in a [jupyter notebook](https://github.com/JanaLasser/agent_based_COVID_SEIRX/blob/dev/school/construct_school_network.ipynb) provided in this repository. Schools are defined by the number of classes they have, the number of students per class, the number of floors these classes are distributed over, and the school type which determines the age structure of the students in the school. A school will have a number of teachers that corresponds to twice the number of classes (which corresponds to approximately the class/teacher ratio in Austrian schools). Every student will have a number of family members drawn from a distribution of household sizes corresponding to Austrian house holds.

In addition to specifying the agent type, nodes also have node attributes that introduce additional parameters into the transmission dynamics: students are part of a ```class``` (```unit```), which largely defines their contact network. Classes are assigned to ```floors``` and have "neighbouring classes" that are situated on the same floor. A small number of random contacts between neighbouring classes are added to the student interaction network, next to the interactions within each class. Teachers have a schedule that specifies the classes they interact with.  

Students within the same class have ```intermediate``` contacts to each other (complete graph). Teachers teaching students in a class have ```intermediate``` contacts to all students in the class. Teachers also have ```intermediate``` contacts to all other teachers (complete graph among teachers), since they regularly meet in the faculty room. Some students have ```far``` contacts to students from other classes, reflecting friendships and passing contacts in hallways to members of other classes. Students also have ```close``` contacts to members of their household. There are no contacts between teachers and family members.  

Students also have an ```age``` that modulates both their ```transmission risk``` and their ```reception risk```. 

## Assumptions
The assumptions made by the model to simplify the dynamics of infection spread and estimates of relevant parameters of virus spread are detailed in the following.

### Parameters
* **Exposure time** (latent time): The time from transmission to becoming infectious is approximated to be $5\pm 1.9$ days ([Ferreti et al. 2020](https://doi.org/10.1126/science.abb6936), [Linton et al. 2020](https://www.mdpi.com/2077-0383/9/2/538), [Lauer et al. 2020](https://www.acpjournals.org/doi/full/10.7326/M20-0504)). Adjust this parameter through the ```exposure_duration``` variable by supplying either a mean duration or the mean and standard deviation of a Weibull distribution.

* **Infectivity duration**: An infected agent is assumed to be infectious for 10.91 $\pm$ 3.95 days after becoming infections ([Walsh et al. 2020](https://doi.org/10.1016/j.jinf.2020.06.067), [You et al. 2020](https://www.sciencedirect.com/science/article/abs/pii/S1438463920302133?via%3Dihub)). Adjust this parameter through the ```infection_duration``` variable by supplying either a mean duration or the mean and standard deviation of a Weibull distribution.

* **Time until symptoms** (incubation time): Humans infected with SARS-CoV2 that develop a clinical course of the disease usually develop symptoms only after they become infectious. We assume the length of the time period between transmission and developing symptoms to be $6.4\pm 0.8$ days ([He et al. 2020](https://www.nature.com/articles/s41591-020-0869-5), [Backer et al. 2020](https://doi.org/10.2807/1560-7917.ES.2020.25.5.2000062)). Adjust this parameter through the ```time_until_symptoms``` variable by supplying either a mean duration or the mean and standard deviation of a Weibull distribution.

* **Infectiousness**: We assume that infectiousness stays constantly high until symptom onset and thereafter decreases monotonically thereafter, until it reaches zero at 10.91 $\pm$ 3.95 days ([He et al. 2020](https://doi.org/10.1038/s41591-020-0869-5), [Walsh et al. 2020](10.1016/j.jinf.2020.06.067)).

* **Symptom probability**: A large proportion of infections with SARS-CoV2 take a subclinical (i.e. asymptomatic) course. We assume that this is true for 40% of infections [(Nikolai et al. 2020)](https://www.sciencedirect.com/science/article/pii/S1201971220307062#bib0100). Nevertheless, a differentiation between residents and personnell might be warranted, with a lower probability to remain asymptomatic for residents, as evidence is mounting that age correlates negatively with the probability to have an asymptomatic course [(McMichael et al. 2020)](https://www.nejm.org/doi/full/10.1056/NEJMoa2005412). Adjust this parameter through the ```symptom_probability``` variable of each agent group.

* **Infectiousness of asymptomatic cases**: We assume that the infectiousness of asymptomatic persons is 40% lower as the infectiousness of symptomatic cases ([Byambasuren et al. 2020](https://www.medrxiv.org/content/10.1101/2020.05.10.20097543v3). Neverthelesss, since their viral load is the same as in symptomatic cases, test sensitivity does not decrease for asymptomatic cases [Walsh et al. 2020](https://doi.org/10.1016/j.jinf.2020.06.067)). Adjust this parameter through the ```subclinical_modifier``` variable.

* **Transmission risk**: For every agent group, a base transmission risk ```transmission_risk``` has to be specified. Transmission risk for a transmission between agents is modulated by the closeness of their contact. The modulation can be specified through the ```infection_risk_contact_type_weights``` – a dictionary, that specifies weights for ```very_far```, ```far```, ```intermediate``` and ```close``` contacts. Transmission risks are also modulated by agent age (if it is specified), to account for lower lung volumes of children [NEEDS A REFERENCE], and by masks (see "mask effectiveness" below). As a note of caution: evidence for the reduced contraction of SARS-CoV-2 in children is still very scarce and the assumptions here are very speculative. Transmission risks in the nursing home model are calibrated such that the outbreak characteristics match the observed empirical data of outbreaks in Austrian nursing homes. In these simulations without interventions, the basic reproduction number R_0 approaches 2.5 to 3 if the index case is introduced through an employee (the value currently reported for SARS-Cov2 spread in the literature, see [Li et al. 2020](https://doi.org/10.1056/NEJMoa2001316), [Wu et al. 2020](http://www.sciencedirect.com/science/article/pii/S0140673620302609) ). If the index case is introduced through a resident, the basic reproduction number lies between 4 and 5, which reflects the confined living conditions and close contacts between residents in nursing homes.

* **Reception risk**: For every agent group, a base reception risk ```reception_risk``` has to be specified. The reception risk is modulated by agent age (if specified), to account for the lower number of receptors responsible for contracting SARS-CoV-2 in children ([Sharif-Askari et al. 2020](https://dx.doi.org/10.1016%2Fj.omtm.2020.05.013)). As a note of caution: evidence for the reduced contraction of SARS-CoV-2 in children is still very scarce and the assumptions here are very speculative.
* **Mask effectiveness**: If agents wear masks, we reduce their transmission risk by 50% and their reception risk by 30%, which are conservative estimates based on the study by [Pan et al. 2020](https://doi.org/10.1101/2020.11.18.20233353).

### Interaction and intervention assumptions
* **Time**: We assume that one model simulation step corresponds to one day. Simulation parameters are chosen accordingly.
* **Tests**: The class ```testing_strategy.py``` implements a variety of different tests, including antigen, PCR and LAMP tests. These tests differ regarding their sensitivity, specificity, the time a test takes until it delivers a result, the time it takes until an infected agent is testable and the time an infected agent stays testable. The test used for diagnostic and preventive testing can be specified at model setup (default is one day turnover PCR test).
* **Quarantine duration**: We assume that agents that were tested positive are isolated (quarantined) for 14 days, according to [recommendations by the WHO](https://www.who.int/publications/i/item/considerations-for-quarantine-of-individuals-in-the-context-of-containment-for-coronavirus-disease-(covid-19)). This time can be changed by supplying the parameter ```quarantine_duration``` to the simulation.
* **Index cases**: There are several ways to introduce index cases to a facility: One way is to introduce a single index case through an agent (specify the agent group through the ```index_case``` parameter) and then simulate the ensuing outbreak in the facility. The second option is to set a probability of an agent to become an index case in each simulation step and choose whether employees, patients or both agent groups can become index cases. To use this option, specify ```index_case='continuous'``` and set the ```index_probability``` parameter to a non-zero risk for the desired agent groups.

## Installation (Linux)
Note: this is currently not working because the dependencies are buggy. This simulation uses the standard scientific python stack (python 3.7, numpy, pandas, matplotlib) plus networkx and mesa. If these libraries are installed, the simulation should work out of the box.

1. Clone the repository:  
```git clone https://github.com/JanaLasser/agent_based_COVID_SEIRX.git```  
Note: if you want to clone the development branch, use  
```git clone --branch dev https://github.com/JanaLasser/agent_based_COVID_SEIRX.git``` 
2. Navigate to the repository  
```cd agent_based_COVID_SEIRX```
3. Create and activate a virtual environment. Make sure you use a Python binary with a version version >= 3.8  
```virtualenv -p=/usr/bin/python3.8 .my_venv```  
```source .my_venv/bin/activate```  
4. Update pip  
``` pip install --upgrade pip```  
5. Install dependencies  
```pip install -r requirements.txt```  

## Running the simulation
The following requires the activation of the virtual environment you created during installation  
```source .my_venv/bin/activate```

I provide exemplary Jupyter Notebooks for [nursing homes](https://github.com/JanaLasser/agent_based_COVID_SEIRX/blob/dev/nursing_home/example_nursing_home.ipynb) and [schools](https://github.com/JanaLasser/agent_based_COVID_SEIRX/blob/dev/school/example_school.ipynb) that illustrate how a simulation model is set up and run for these two applications, how results are visualised and how data from a model run can be collected. Run the example notebook from the terminal:  
```jupyter-notebook nursing_home/example_nursing_home.ipynb```  

or  

```jupyter-notebook school/example_school.ipynb```

I also provide the [Jupyter Notebook](https://github.com/JanaLasser/agent_based_COVID_SEIRX/blob/dev/nursing_home/screening_frequency_data_creation.ipynb) used to run the simulations and create the data used in the publication **Agent-based simulations for optimized prevention of the spread of SARS-CoV-2 in nursing homes** as well as the [Jupyter Notebook](https://github.com/JanaLasser/agent_based_COVID_SEIRX/blob/master/nursing_home/screening_frequency_analysis.ipynb) used to create the heatmaps for the different analysed screnarios from the simulation data. Run these notebooks from the terminal using:

```jupyter-notebook nursing_home/screening_frequency_data_creation.ipynb```  

and  

```jupyter-notebook nursing_home/screening_frequency_analysis.ipynb```  


## Acknowledgements
I would like to thank [Peter Klimek](https://www.csh.ac.at/researcher/peter-klimek/) from Complexity Science Hub Vienna and Thomas Wochele-Thoma from [Caritas Austria](https://www.caritas.at/) for the fruitful discussions that led to the development of this project.

## References
Ferretti, Luca, et al. "Quantifying SARS-CoV-2 transmission suggests epidemic control with digital contact tracing." Science 368.6491 (2020). [DOI: 10.1126/science.abb6936](https://doi.org/10.1126/science.abb6936)

Linton, N. M., Kobayashi, T., Yang, Y., Hayashi, K., Akhmetzhanov, A. R., Jung, S. M., ... & Nishiura, H. (2020). Incubation period and other epidemiological characteristics of 2019 novel coronavirus infections with right truncation: a statistical analysis of publicly available case data. Journal of clinical medicine, 9(2), 538. [DOI: 10.3390/jcm9020538](https://doi.org/10.3390/jcm9020538)  

Lauer, S. A., Grantz, K. H., Bi, Q., Jones, F. K., Zheng, Q., Meredith, H. R., ... & Lessler, J. (2020). The incubation period of coronavirus disease 2019 (COVID-19) from publicly reported confirmed cases: estimation and application. Annals of internal medicine, 172(9), 577-582. [DOI: 10.7326/M20-0504](https://doi.org/10.7326/M20-0504)  

Walsh, K. A., Jordan, K., Clyne, B., Rohde, D., Drummond, L., Byrne, P., ... & O'Neill, M. (2020). SARS-CoV-2 detection, viral load and infectivity over the course of an infection: SARS-CoV-2 detection, viral load and infectivity. Journal of Infection. [DOI: 10.1016/j.jinf.2020.06.067](10.1016/j.jinf.2020.06.067)  

You, C., Deng, Y., Hu, W., Sun, J., Lin, Q., Zhou, F., ... & Zhou, X. H. (2020). Estimation of the time-varying reproduction number of COVID-19 outbreak in China. International Journal of Hygiene and Environmental Health, 113555. [DOI: 10.1016/j.ijheh.2020.113555](https://doi.org/10.1016/j.ijheh.2020.113555)

Backer Jantien A, Klinkenberg Don, Wallinga Jacco. Incubation period of 2019 novel coronavirus (2019-nCoV) infections among travellers from Wuhan, China, 20–28 January 2020. Euro Surveill. 2020;25(5):pii=2000062. [DOI: 10.2807/1560-7917.ES.2020.25.5.2000062](https://doi.org/10.2807/1560-7917.ES.2020.25.5.2000062)

He, X., Lau, E. H., Wu, P., Deng, X., Wang, J., Hao, X., ... & Mo, X. (2020). Temporal dynamics in viral shedding and transmissibility of COVID-19. Nature medicine, 26(5), 672-675. [DOI: 10.1038/s41591-020-0869-5](https://doi.org/10.1038/s41591-020-0869-5)  

Nikolai, L. A., Meyer, C. G., Kremsner, P. G., & Velavan, T. P. (2020). Asymptomatic SARS Coronavirus 2 infection: Invisible yet invincible. International Journal of Infectious Diseases. [DOI: 10.1016/j.ijid.2020.08.076](https://doi.org/10.1016/j.ijid.2020.08.076)  

McMichael, T. M., Currie, D. W., Clark, S., Pogosjans, S., Kay, M., Schwartz, N. G., ... & Ferro, J. (2020). Epidemiology of Covid-19 in a long-term care facility in King County, Washington. New England Journal of Medicine, 382(21), 2005-2011. [DOI: 10.1056/NEJMoa2005412](https://doi.org/10.1056/NEJMoa2005412)  

Li, Q., Guan, X., Wu, P., Wang, X., Zhou, L., Tong, Y., ... & Xing, X. (2020). Early transmission dynamics in Wuhan, China, of novel coronavirus–infected pneumonia. New England Journal of Medicine. [DOI: 10.1056/NEJMoa2001316](https://doi.org/10.1056/NEJMoa2001316)  

Wu, J. T., Leung, K., & Leung, G. M. (2020). Nowcasting and forecasting the potential domestic and international spread of the 2019-nCoV outbreak originating in Wuhan, China: a modelling study. The Lancet, 395(10225), 689-697. [DOI: 10.1016/S0140-6736(20)30260-9](https://doi.org/10.1016/S0140-6736(20)30260-9)  

Pan, J., Harb, C., Leng, W., & Marr, L. C. (2020). Inward and outward effectiveness of cloth masks, a surgical mask, and a face shield. medRxiv. [DOI: 10.1101/2020.11.18.20233353](https://doi.org/10.1101/2020.11.18.20233353)

Sharif-Askari, N. S., Sharif-Askari, F. S., Alabed, M., Temsah, M. H., Al Heialy, S., Hamid, Q., & Halwani, R. (2020). Airways Expression of SARS-CoV-2 Receptor, ACE2, and TMPRSS2 Is Lower in Children Than Adults and Increases with Smoking and COPD. Mol Ther Methods Clin Dev, 18, 1-6. [DOI: 10.1016%2Fj.omtm.2020.05.013](https://dx.doi.org/10.1016%2Fj.omtm.2020.05.013)
