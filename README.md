# Agent based simulation of the spread of COVID-19 in confined spaces
**Author: Jana Lasser, Complexity Science Hub Vienna (lasser@csh.ac.at)**

A simulation tool to explore the spread of COVID-19 in small communities such as nursing homes or schools via agent-based modeling (ABM) and the impact of prevention measures. The model follows an SEIRX approach, building on the agent based simulation framework [mesa](https://mesa.readthedocs.io/en/master/) in which agents can be susceptible (S), exposed (E), infected (I), removed (R) or quarantined (X) and is based on explicitly defined and dynamic contact networks between agents. The model offers the possibility to explore the effectiveness of various testing, tracing and quarantine strategies and other interventions such as ventilation and mask-wearing.

![](https://github.com/JanaLasser/agent_based_COVID_SEIRX/blob/packaging/src/scseirx/img/nursing_home_contact_network_illustration.png?raw=true)

*This software is under development and intended to respond rapidly to the current situation. Please use it with caution and bear in mind that there might be bugs.*


Reference:  

_Lasser, J. (2020). Agent based simulation of the spread of COVID-19 in nursing homes. [DOI](https://doi.org/10.5281/zenodo.4613202): 10.5281/zenodo.4613202_

## Table of contents
* [Simulation Design](#simulation-design)
    * [Actors](#actors)
    * [Transmissions](#transmissions)
    * [Containment strategies](#containment-strategies)
* [Implementation](#implementation)
    * [Model](#model)
    * [Agents](#agents)
    * [Testing](#testing)
    * [Additional modules](#additional-modules)
* [Applications](#applications)
    * [Nursing Homes](#nursing-homes)
    * [schools](#schools)
* [Assumptions and Approximations](#assumptions-and-approximations)
    * [Epidemiological parameters](#epidemiological-parameters)
    * [Intervention Measure Effectiveness](#intervention-measure-effectiveness)
* [Calibration](#calibration)
    * [Household Transmissions](#household-transmissions)
    * [Calibration for Schools](#calibration-for-schools)
    * [Calibration for Nursing Homes](#calibration-for-nursing-homes)
* [Installation Linux](#installation-linux)
* [Running the simulation](#running-the-simulation)
* [Acknowledgements](#acknowledgements)
* [References](#references)

## Simulation Design

### Actors
#### States
Simulations can have several types of actors (agents), for example residents & employees in nursing homes, or teachers, students and family members in schools. Infections are introduced through agents that have a certain probability to become an index case or that are explicitly chosen as index case in the beginning of the simulation. The contact network defines which agents interact with which other agents on which day of the week. Different contact types modulate the infection transmission risk between "close" and "very far". Contact networks are stored as a [networkx](https://networkx.org/) graph with edge attributes for different contact types. In every step (day) of the simulation, agents interact according to the contact network and can transmit the infection to the agents they are in contact with. Depending on its infection state, an agent has one of four states: susceptible (S), exposed (E), infected (I) or removed (R). In addition, agents can be quarantined (X) and develop symptoms when they are infected. While infected (I), agents can by symptomatic or asymptomatic. 

#### Attributes
Next to the dynamic states pertaining to the infection state of an agent, agents also have (static) attributes that influence how they interact with the infection dynamics:
* Exposure duration: The time between transmission and becoming infectious (exposure duration) for every agent is drawn from a Weibull distribution and might be different for every agent (see section [Epidemiological Parameters](#epidemiological-parameters)) for details).
* Time until symptoms: The time between transmission and (potentially) showing symptoms for every agent is drawn from a Weibull distribution and might be different for every agent (see section [Epidemiological Parameters](#epidemiological-parameters) for details).
* Infection duration: The time between transmission and ceasing to be infectious is drawn from a Weibull distribution and might be different for every agent (see section [Epidemiological Parameters](#epidemiological-parameters) for details).
* Age: especially in the school setting, where children are involved, the age of the agents plays an important role, since children have a somewhat reduced risk to transmit or receive an infection (see section [Transmissions](#transmissions) for details).

#### Viral Load
An infection with SARS-CoV-2 causes an infected person's body to replicate the virus and the viral load in the person will vary, depending on the progression of the infection. In this simulation, viral load influences two distinct processes: the ability of testing technologies to detect an infection (detection threshold, see section [Test Technologies](#test-technologies)) and the infectiousness of agents (see section [Transmissions](#transmissions)). So far, we do not model the dependence of viral load on time explicitly for every agent. As pertaining to infectiousness, we model the viral load as trapezoid function that is high at the beginning and drops to zero over the course of the infection. For different testing technologies, we model the time between transmission and detection threshold individually for every testing technology used.

### Transmissions
Transmissions are modeled as [Bernoulli-Trials](https://en.wikipedia.org/wiki/Bernoulli_trial) with a probability p of success and a probability q = 1 - p of failure. In every step of the simulation (one step corresponds to one day), every infected and non-quarantined agent performs this Bernoulli trial once for every other (non-infected and non-quarantined) agent they are in contact with. The overall probability of a successful transmission between two agents is reduced by a range of factors that reflect both biological factors and intervention measures that can act on both the _transmitter_ of the infection as well as the _receiver_ of the infection. If q_m is the probability of failure of transmission due to wearing a mask and b is the baseline transmission risk without any modifications, then the modified probability of successful transmission between two agents is

p = 1 - [1 - b(1 - q_m)].

We currently account for eight different factors that can influence the transmission risk in different settings (see sections [Epidemiological Parameters](#epidemiological-parameters) and [Intervention measure effectiveness](#intervention-measure-effectiveness) for details on how these factors can be set and how values for them are chosen):
* q_1: Modification of the transmission risk due to the type of contact between agents. Here, q_1 = 1 for a household contact (contact type "close", no modification), whereas q_1 < 1 for other contact types ("intermediate", "far", "very far", reduction of transmission risk). The value of q_1 depending on the type of the contact has to be specified via the model parameter ```infection_risk_contact_type_weights``` at model setup. Contact types between each two agents are stored in the contact network supplied to the model.
* q_2: Modification of the transmission risk due to the age of the _transmitting agent_. The dependence of transmission risk on age is set via the model parameter ```age_symptom_discount``` at model setup.
* q_3: Modification of the reception risk due to the age of the _receiving agent_. The dependence of reception risk on age is approximated to be the same as the dependence of the transmission risk on age and is therefore also set via the model parameter ```age_symptom_discount``` at model setup.
* q_4: Modification of the transmission risk due to the progression of the infection. This dependence is currently hard-coded, based on literature values.
* q_5: Modification of the transmission risk due to the type of the course of the infection (symptomatic, asymptomatic). This parameter is set via the model parameter ```subclinical_modifier``` at model setup.
* q_6: Modification of the transmission risk due to mask wearing of the _transmitting agent_. This parameter is set via the model parameter ```mask_filter_efficiency["exhale"]``` at model setup. Whether or not an agent group is wearing masks has to be specified via the ```agent_types["mask"]``` parameter and the contact types that are affected by mask-wearing are hard-coded in the model (for example, household contacts are generally not affected by mask-wearing).
* q_7: Modification of the reception risk due to mask wearing of the _receiving agent_. This parameter is set via the model parameter ```mask_filter_efficiency["inhale"]``` at model setup. Whether or not an agent group is wearing masks has to be specified via the ```agent_types["mask"]``` parameter and the contact types that are affected by mask-wearing are hard-coded in the model (for example, household contacts are generally not affected by mask-wearing).
* q_8: Modification of the transmission risk due to room ventilation. This parameter is set via the model parameter ```transmission_risk_ventilation_modifier``` at model setup.  
* q_9: Modification of the transmission risk due to vaccination. This parameter is set via the midel parameter ```transmission_risk_vaccination_modifier``` at model setup.

The baseline transmission risk is set via the model parameter ```base_transmission_risk``` at model setup.

Therefore, for example in a school setting where agents are wearing masks, rooms are ventilated and the age of agents is important for the transmission dynamics, the overall success probability for a transmission is defined as  

p = 1 - [1 - b(1 - q_1)(1 - q_2)(1-q_3)(1 - q_4)(1 - q_5)(1 - q_6)(1 - q_7)(1 - q_8)(1 - q_9)].

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
The simulation consists of a _model_ that stores model parameters, the agent contact network and references to all agents. In every step of the simulation, the model initiates agent interactions and executes the testing and tracing strategy as well as data collection. Model parameters and parameters for the testing strategy, as well as a specification of the agent types and their respective parameters have to be passed to the  model instance at time of creation, if values other than the specified default values (see sections [Epidemiological Parameters](#epidemiological-parameters) and [Intervention measure effectiveness](#intervention-measure-effectiveness)) should be used. The model is implemented as a class (```model_SEIRX```) that inherits from [mesa's](https://mesa.readthedocs.io/en/stable/) ```Model``` and implements a ```step()``` function that is called in every time-step. Every scenario (so far: nursing homes and schools) implements its own model class which inherits from ```model_SEIRX```, where functionality deviating from the behaviour specified in ```model_SEIRX``` is implemented. This can for example a custom ```step()``` or ```calculate_transmission_probability()``` function and specify scenario-specific data collection functions.

Test technologies for preventive screening and diagnostic testing (diagnostic tests, background screens and follow-up screens) can be specified separately, using the model parameters ```diagnostic_test_type``` and ```preventive_screening_test_type```.

### Agents
Similarly to the model, agents have a base-class ```agent_SEIRX``` that inherits from [mesa's](https://mesa.readthedocs.io/en/stable/) ```Agent```. The agent class implements agent states and counters and functions necessary for simulating contacts between agents and advancing states. Different agent types needed in the scenarios inherit from this base class and might implement additional functionality. Currently there are five agent types: resident and employee (nursing home scenario) and teacher, student and family member (school scenario). These agent types are implemented in separate classes which inherit from the agent base-class. 

### Testing
The testing strategy is contained in ```testing_strategy.py```, a class different from the SEIRX base model but is created with parameters passed through the SEIRX constructor. This is to keep parameters and information related to testing and tracing in one place, separate from the infection dynamics model. The testing class also stores information on the sensitivity, specificity and turnover time of a range of tests and can be easily extended to include additional testing technologies.

### Additional modules
* The module ```analysis_functions.py``` provides a range of functions to analyse data from model runs.
* The module ```viz.py``` provides some custom visualization utility to plot infection time-lines and agent states on a network, given a model instance.
* the module ```school/construct_school_network.py``` provides functionality to construct contact networks for different school types

## Applications
### Nursing homes
Nursing homes implement agents types ```resident``` and ```employee```, as well as the ```model_nursing_home``` (all located in the ```nursing_home``` sub-folder).  

The contact networks for nursing homes are specified through resident relations in the homes (room neighbors, table neighbors at joint meals, other shared areas). Employees interact with all other employees and all residents in the same living quarter of the nursing home every day. We provide several exemplary contact networks (see ```data/nursing_homes```), representing different architectures of nursing homes with different numbers of living quarters. These contact networks are abstractions of empirically determined interaction relations in Austrian nursing homes. By default, resident roommates are defined as "close contact", residents that share tables during joint meals as well as all resident - employee contacts are considered "intermediate", and contacts between residents that only share the same living quarters but not the same room or table are considered to be "far". See the [example notebook](https://github.com/JanaLasser/agent_based_COVID_SEIRX/blob/master/example_nursing_home.ipynb) for an exemplary simulation of a nursing home scenario.

### Schools
Schools implement agent types ```teachers```, ```students``` and ```family_members``` of students, as well as the ```model_school``` (all located in the ```school``` sub-folder). See the [example notebook](https://github.com/JanaLasser/agent_based_COVID_SEIRX/blob/master/example_school.ipynb) for an exemplary simulation of a school scenario.

The contact networks for schools are modeled to reflect common structures in Austrian schools (see the [school-type specific documentation](https://github.com/JanaLasser/agent_based_COVID_SEIRX/blob/master/school/school_type_documentation.ipynb) for details). Schools are defined by the number of classes they have, the number of students per class, the number of floors these classes are distributed over, and the school type which determines the age structure of the students in the school. A school will have a number of teachers that is determined automatically, depending on the number of classes and the school type. Every student will have a number of family members drawn from a distribution of household sizes corresponding to Austrian house holds. For this application, we construct three distinct types of network:
* **Representative schools**: We construct networks for the seven most common school types in Austria: Primary schools, primary schools with daycare, lower secondary schools, lower secondary schools with daycare, upper secondary schools, secondary schools (lower & upper secondary education) and secondary schools with daycare (only for students in the lower secondary age bracket). These schools are constructed with characteristics (number of classes, number of students per class) representing the average over all the existing schools of a given school type in Austria (statistics from the year 2017/18). See the respective [Jupyter notebook](https://github.com/JanaLasser/agent_based_COVID_SEIRX/blob/master/school/construct_school_networks_representative_schools.ipynb) for details.
* **Calibration schools**: These contact networks are similar to those of the representative schools described above. There are only two differences: Firstly, we remove all household members that are not siblings that go to the same school from the households, since we do not need them for the calibration, and the simulations run faster on smaller networks. Secondly, we create 500 instances of each contact network for each school type. Sibling connections are generated randomly, based on the households that are generated for each child. To sample a range of sibling constellations that enable contacts between different classes in the same school, we create a number of networks equal to the size of the ensembles simulated for each parameter combination during calibration. See the respective [Jupyter notebook](https://github.com/JanaLasser/agent_based_COVID_SEIRX/blob/master/school/construct_school_network_calibration_schools.ipynb) for details.
* **Various school layouts**: In the interactive visualization of the infection dynamics, we enable the user to explore schools with different characteristics (number of classes, number of students per class). We create the contact networks for all the different possible combinations of school characteristics. See the respective [Jupyter notebook](https://github.com/JanaLasser/agent_based_COVID_SEIRX/blob/master/school/construct_school_network_all_layouts.ipynb) for details.


## Assumptions and Approximations
The assumptions and approximations made by the model to simplify the dynamics of infection spread and estimates of relevant epidemiological parameters are detailed in the following. 

* **Time**: We assume that one model simulation step corresponds to one day. Simulation parameters are chosen accordingly.

* **Index cases**: There are several ways to introduce index cases to a facility: One way is to introduce a single index case through an agent (specify the agent group through the ```index_case``` parameter) and then simulate the ensuing outbreak in the facility. The agent will be randomly chosen from all agents in the corresponding group. The second option is to set a probability of an agent to become an index case in each simulation step and choose whether employees, patients or both agent groups can become index cases. To use this option, specify ```index_case='continuous'``` and set the ```index_probability``` in the ```agent_types``` parameter to a non-zero risk for the desired agent groups.

### Epidemiological parameters
* **Exposure duration** (latent time): The time from transmission to becoming infectious is approximated to be 5 +- 1.9 days ([Ferreti et al. 2020](https://doi.org/10.1126/science.abb6936), [Linton et al. 2020](https://www.mdpi.com/2077-0383/9/2/538), [Lauer et al. 2020](https://www.acpjournals.org/doi/full/10.7326/M20-0504)), which we use in both scenarios and as default value for the simulation. Adjust this parameter through the ```exposure_duration``` variable by supplying either a mean duration or the mean and standard deviation of a Weibull distribution. If mean & standard deviation of a Weibull distribution are supplied, the exposure duration of every agent is drawn from this distribution, such that 0 < exposure duration <= time until symptoms and exposure duration < infection duration.

* **Infection duration**: An infected agent is assumed to be infectious for 10.91 +- 3.95 days after becoming infections ([Walsh et al. 2020](https://doi.org/10.1016/j.jinf.2020.06.067), [You et al. 2020](https://www.sciencedirect.com/science/article/abs/pii/S1438463920302133?via%3Dihub)), which we use in both scenarios and as default value for the simulation. Adjust this parameter through the ```infection_duration``` variable by supplying either a mean duration or the mean and standard deviation of a Weibull distribution. If mean & standard deviation of a Weibull distribution are supplied, the infection duration of every agent is drawn from this distribution, such that infection duration > exposure duration.

* **Time until symptoms** (incubation time): Humans infected with SARS-CoV2 that develop a clinical course of the disease usually develop symptoms only after they become infectious. We assume the length of the time period between transmission and developing symptoms to be 6.4 +- 0.8 days ([He et al. 2020](https://www.nature.com/articles/s41591-020-0869-5), [Backer et al. 2020](https://doi.org/10.2807/1560-7917.ES.2020.25.5.2000062)), which we use in both scenarios and as default value for the simulation. Adjust this parameter through the ```time_until_symptoms``` variable by supplying either a mean duration or the mean and standard deviation of a Weibull distribution. If mean & standard deviation of a Weibull distribution are supplied, the time until symptoms of every agent is drawn from this distribution, such that time until symptoms >= exposure duration.

* **Infectiousness**: We assume that infectiousness stays constantly high until symptom onset and thereafter decreases monotonically thereafter, until it reaches zero at the end of the infection duration ([He et al. 2020](https://doi.org/10.1038/s41591-020-0869-5), [Walsh et al. 2020](10.1016/j.jinf.2020.06.067)). We model this by a trapezoid function that is one (i.e. infectiousness is equal to the base transmission risk, if no other modifiers apply) in the time between the end of the exposure duration and the onset of symptoms. After symptom onset, infectiousness declines linearly over a number of days equal to infection duration - time until symptoms. This behaviour is currently hard-coded (see ```get_transmission_risk_progression_modifier()``` in the ```model_SEIRX``` class).

![](https://github.com/JanaLasser/agent_based_COVID_SEIRX/blob/packaging/src/scseirx/img/transmission_risk_progression.png?raw=true)

* **Symptom probability**: A large proportion of infections with SARS-CoV2 take a subclinical (i.e. asymptomatic) course. For the school scenario, we use the empirically observed relation (from cluster tracing by AGES in Austria) between the number of symptomatic courses and the age of the children: at age 18 or above, the probability to develop a symptomatic course is approximately 85%. For younger children, this probability drops by about 2.5% for every year the child is younger than 18. This information is supplied to the model through the parameter ```age_symptom_discount```, which defines a line where the symptom probability for adults is the intercept and the slope is the amount by which this probability is decreased for every year below 18. For the nursing home scenario, we assume that about 40% of infections remain asymptomatic, since this is the information found in the literature that is best matching our scenario ([(Nikolai et al. 2020)](https://www.sciencedirect.com/science/article/pii/S1201971220307062#bib0100)), and we do not have our own empirical data for this scenario. Nevertheless, a differentiation between residents and personnel might be warranted, with a lower probability to remain asymptomatic for residents, as evidence is mounting that age correlates negatively with the probability to have an asymptomatic course in this age group as well [(McMichael et al. 2020)](https://www.nejm.org/doi/full/10.1056/NEJMoa2005412).

![](https://github.com/JanaLasser/agent_based_COVID_SEIRX/blob/packaging/src/scseirx/img/age_symptom_discount.png?raw=true)

* **Infectiousness of asymptomatic cases**: We assume that the infectiousness of asymptomatic persons is 40% lower than the infectiousness of symptomatic cases ([Byambasuren et al. 2020](https://www.medrxiv.org/content/10.1101/2020.05.10.20097543v3). Neverthelesss, since their viral load is the same as in symptomatic cases, test sensitivity does not decrease for asymptomatic cases ([Walsh et al. 2020](https://doi.org/10.1016/j.jinf.2020.06.067)). Adjust this parameter through the ```subclinical_modifier``` variable (default=0.6).

* **Base transmission risk**: The base transmission risk specifies the risk of transmitting an infection during a contact of type "close" (considered a household contact). We choose the base transmission risk in such a way that household members have a risk of 37.8% of contracting the disease if another household member is infected, reflecting the current state of the literature on household transmissions of SARS-CoV-2 ([Madwell et al. 2020](https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2774102)). This yields a base transmission risk of 0.74 or about 7.4% per day between household contacts.

* **Reception risk age modification**: There is mounting evidence that the risk of a susceptible agent to "receive" an infection from an infected agent depends on age, since younger children seem to express fewer receptors responsible for admitting the virus to the body [Sharif-Askari et al. 2020](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7242205/). We model this together with a reduced risk of transmission (see below) for younger people. This is controlled via the parameter ```age_transmission_risk_discount```, which defines a ```slope``` that specifies the reduction in reception risk per year an agent is younger than 18. This parameter is calibrated to match the empirical observations of infection dynamics.

* **Transmission risk age modification**: Transmission risks are also modulated by agent age, to reflect our empirical observations that the infection dynamics seem to be weaker in smaller children. This might be due to the higher number of asymptomatic cases among children. We would therefore expect less coughing and, together with the smaller lung volumes of children, an overall weaker expulsion of aerosols. Reliable literature on this topic is still missing, therefore we decided to only use one parameter for both the reception risk age modification and the transmission risk age modification (```age_transmission_discount```), which is calibrated for different settings.

![](https://github.com/JanaLasser/agent_based_COVID_SEIRX/blob/packaging/src/scseirx/img/age_transmission_discount.png?raw=true)

* **Infection risk contact type modification**: Transmission risk for a transmission between agents is modulated by the closeness of their contact. The modulation can be specified through the parameter ```infection_risk_contact_type_weights``` – a dictionary, that specifies weights for ```very_far```, ```far```, ```intermediate``` and ```close``` contacts. By definition, the weight for ```close``` should be set to one, since this contact type corresponds to household contacts. For our application scenarios, these weights are calibrated such that they reflect the empirically observed infection dynamics in a given scenario. Weights in the nursing home model are calibrated such that the outbreak characteristics match the observed empirical data of outbreaks in Austrian nursing homes. In these simulations without interventions, the basic reproduction number R_0 approaches 2.5 to 3 if the index case is introduced through an employee (the value currently reported for SARS-Cov2 spread in the literature, see [Li et al. 2020](https://doi.org/10.1056/NEJMoa2001316), [Wu et al. 2020](http://www.sciencedirect.com/science/article/pii/S0140673620302609) ). If the index case is introduced through a resident, the basic reproduction number lies between 4 and 5, which reflects the confined living conditions and close contacts between residents in nursing homes.


### Intervention measure effectiveness
* **Mask effectiveness**: If agents wear masks, we reduce their transmission risk by 50% and their reception risk by 30%, which are conservative estimates for the effectiveness of surgical masks based on the study by [Pan et al. 2020](https://doi.org/10.1101/2020.11.18.20233353). These mask filter efficiencies for transmission and reception (```exhale``` and ```inhale```) can be modified through the parameter ```mask_filter_efficiency```.

* **Ventilation**: For the school scenario, we calculate the ventilation efficiency of short and intensive ventilation of the classroom once or twice per hour during one lesson for using the [ventilation efficiency calculator](https://www.mpic.de/4747361/risk-calculator). Room areas for classrooms were assumed to follow the [Bauverordnung §5](https://www.ris.bka.gv.at/GeltendeFassung.wxe?Abfrage=LrBgld&Gesetzesnummer=10000209): min. 1.6 m² / student and min 50² in primary and lower secondary schools. Since the maximum class size we simulate (secondary schools) is 30 students (30*1.6 m² = 48m²), we can assume all classrooms have a size of approximately 50 m². The individual infection risk is independent of the number of people in the room. Ventilation efficiencies were calculated for both a teacher and student index case and for mask-wearing of transmitting and receiving agents, but ventilation efficiency does not depend on either of these parameters. We calculated transmission risk reductions for an air exchange rate of 2 (corresponding to one short and intense ventilation per hour) and an air exchange rate of 4 (corresponding to two short and intense ventilations per hour). The corresponding reductions in transmission risk are 64% (one ventilation) and 80% (2 ventilations) respectively. Ventilation efficiencies can be adjusted through the ```transmission_risk_ventilation_modifier``` parameter.

* **Tests**: The class ```testing_strategy.py``` implements a variety of different tests, including antigen, PCR and LAMP tests. These tests differ regarding their sensitivity, specificity, the time a test takes until it delivers a result, the time it takes until an infected agent is testable and the time an infected agent stays testable. The test used for diagnostic and preventive testing can be specified at model setup (default is one day turnover PCR test).

* **Quarantine duration**: We assume that agents that were tested positive are isolated (quarantined) for 14 days, according to [recommendations by the WHO](https://www.who.int/publications/i/item/considerations-for-quarantine-of-individuals-in-the-context-of-containment-for-coronavirus-disease-(covid-19)). This time can be changed by supplying the parameter ```quarantine_duration``` to the simulation.

* **Vaccination**: Using the the inout dictionary ```vaccination agent``` the probabilities for the different agent types can be modified. A default of 10 percent for every agent type is assumed. The ```transmission_risk_vaccination_modifier``` is set to 1 by default which means that vaccination has no impact on the transmission of an infected agent to a susceptible, vaccinated agent. It is recommended to set the ```transmission_risk_vaccination_modifier``` to 0.95 for assessing the effects of the vaccines BNT162b2 (BioN-Tech/Pfizer) and mRNA-1273 (Moderna) [Polack et al. 2020](https://www.nejm.org/doi/10.1056/NEJMoa2034577), [Voysey et al. 2020](https://www.thelancet.com/action/showPdf?pii=S0140-6736%2820%2932661-1).

## Calibration
The most important part of any agent based model is its calibration. As described above, the model has many parameters that can be set and will influence the dynamics of infection spread. Some parameter choices can be based on existing literature (such as the effectiveness of masks or ventilation) or directly observable characteristics of infection spread in our settings (such as the age dependence of the probability to develop a symptomatic course). Depending on the setting, there will be a number of free parameters in the model that have to be calibrated to reproduce the observed dynamics of infection spread as closely as possible. In our application, these are a total of three parameters for nursing homes and four parameters for schools:
1. **The base transmission risk of a household contact**, i.e. the probability to transmit an infection through a contact of type "close". This parameter is calibrated for both the nursing home and school setting.
2. **The weight of "intermediate" contacts**, i.e. how much the transmission risk through an intermediate contact is reduced as compared to a close contact. If a close contact has a weight of 1, a weight of 0.8 for an intermediate contact means, that the transmission risk through an intermediate contact is reduced by 20%. This parameter is calibrated for both the nursing home and school setting.
3. **The weight of "far" contacts**, i.e. how much the transmission risk through a far contact is reduced as compared to a close contact. If a close contact has a weight of 1, a weight of 0.4 for an intermediate contact means, that the transmission risk through an intermediate contact is reduced by 60%. This parameter is calibrated for both the nursing home and school setting.
4. **The reduction of transmission and reception risk for children**. Here we have to calibrate a linear discount factor for every year a child is younger than 18. If the discount factor is 0.02, then children aged 6 will have a transmission risk that is reduced by 24% as compared to people aged 18 and above. This parameter is only calibrated for the school setting, since this is the only setting were children are involved, for which significantly reduced transmission and reception risks are assumed.

### Household transmissions
Most recent research indicates that there is a 37.8% risk of adult members of the same household than an infected person to get infected themselves over the course of the infection ([Madwell et al. 2020](https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2774102)). Our aim is to calibrate the base transmission risk _b_ between adult agents in our systems to reflect this empirical finding. Transmission risk is influenced by a range of factors (see section [Transmissions](#transmissions)). For household transmissions between adults, the only relevant factors that influence the transmission risk is the reduction due to progression of the disease _q_4(t)_ (in later stages of the infection the transmission risk decreases) and the reduction in case of an asymptomatic course _q_5_. The values for both of these factors are taken from literature. For one contact on day _i_ the probability of a successful transmission (transmission risk) is therefore

_p_i = 1 - [1 - b(1 - q_4_i)(1 - q_5)]_

In our model, we draw the relevant epidemiological parameters (exposure duration, infection duration, symptomatic course) individually for every agent from distributions. To calibrate, we create pairs of agents and let one of them be infected. We then simulate the whole course of the infection (from day 0 to the end of the infection duration) and check for a transmission with probability of success _p_i_ on every day _i_. We minimize the difference between the expected number of successful infections (37.8%) and the simulated number of successful infections by varying _b_. This results in a value of _b=0.074_ or an average risk of 7.4% for a household member per day to get infected.  
We note that the reduction of transmission risk and reception risk due to age of the transmitting and receiving agents is treated separately. This is why we only calibrate the transmission risk between adults here and calibrate the age discount factor separately.

### Calibration for schools
For schools, we simultaneously calibrate the following parameters
* weight of the contact type "intermediate" (```infection_risk_contact_type_weights['intermediate']```),
* weight of the contact type "far" (```infection_risk_contact_type_weights['far']```) and
* the age transmission risk discount (```age_transmission_risk_discount['slope']```).

For the calibration, we use observations of SARS-CoV-2 outbreaks in Austrian schools in the weeks 35-46 of the year 2020. We optimize two distinct target observables:
1. The distribution of outbreak sizes
2. The distribution of cases to the agent groups "teacher" and "student".
For the optimization, we calculate the sum of the [Chi-squared distances](https://link.springer.com/referenceworkentry/10.1007%2F978-0-387-32833-1_53) between the empirically observed distributions and the distributions generated from ensembles of 500 runs for every parameter combination. 

For the other simulation parameters, we use settings that most closely match the situation in Austrian schools in the period of time from which the empirical observations stem: 
* Index cases are drawn from the empirically observed distribution of index cases between teachers and students.
* The age dependence of the probability of developing a symptomatic course is matched to the empirically observed age dependence.
* Only diagnostic testing with PCR tests with a one-day turnover was in place, followed by a background screen in case of a positive result.
* There were no preventive screens and no follow-up tests after a background screen.
* Contacts of type "close" and "intermediate" were considered to be K1 contacts and were quarantined for 10 days in case of a positive test result and there was no "liberating testing" in place.
* Teachers and students did not regularly wear masks during lessons.
* Teachers and students did wear masks in hallways and shared community areas and contacts between students of different classes were avoided.
* All students of a class were present every day (i.e. no halving of classes).

Using these settings, to find optimal values for the parameters, we first conduct a random search in the parameter grid spanned by the following ranges:
* contact weight intermediate: [0:1:0.05],
* contact weight far: [0:1:0.05] and
* age transmission discount: [-0.1:0:0.02],
where we impose the additional constraint on parameter combinations that contact weight intermediate > contact weight far. We randomly choose 100 parameter combinations out of the 950 possible parameter combinations and simulate ensembles of 500 runs for each parameter combination and school type. For school characteristics we choose a number of classes and students per class that most closely matches the average characteristics of an Austrian school of the given school type (see section [Schools](#schools)). Since the empirical data we have does not differentiate between schools with and without daycare of a given school type, we assume that 50% of the schools of a given school type are schools with daycare. This approximates the [percentage of schools with daycare in Austria](https://www.kdz.eu/de/content/fact-sheets-pflichtschule-und-tagesbetreuung). The [Austrian school statistics](https://www.bmbwf.gv.at/Themen/schule/schulsystem/gd.html) also do not differentiate between schools with and without daycare. Therefore we assume that schools with and without daycare are not significantly different in the number of classes and number of students per class. We therefore simulate ensembles for primary schools, primary schools with daycare, lower secondary schools, lower secondary schools with daycare, upper secondary schools (no daycare in this school type), secondary schools and secondary schools with daycare with the following characteristics:

school type     | # classes | # students / class
--------------- | --------- | ------------------
primary         | 8         | 19
lower secondary | 8         | 18
upper secondary | 10        | 23
secondary       | 28        | 24

We calculate the overall difference between the simulated distribution of outbreak sizes and infected agent types as sum of the Chi-squared distances of every school type, weighted by the number of empirical observations for this school type.

After we identify the parameter combination that minimizes this Chi-squared distance in the random grid search, we perform a refined grid search around the current optimal parameter combination and repeat the optimization process as described above. We find that a parameter combination of
* contact weight intermediate = 0.85
* contact weight far = 0.75
* age transmission discount = -0.02
produces outbreak characteristics that most closely match the empirically observed outbreaks. We use these parameter combination for all subsequent simulations to analyze the effect of different prevention strategies on outbreak characteristics in schools.

### Calibration for nursing homes
TODO


## Installation Linux
1. Clone the repository:  
```git clone https://github.com/JanaLasser/agent_based_COVID_SEIRX.git```  
Note: if you want to clone the development branch, use  
```git clone --branch dev https://github.com/JanaLasser/agent_based_COVID_SEIRX.git``` 
2. Navigate to the repository  
```cd agent_based_COVID_SEIRX```
3. Create and activate a virtual environment. Make sure you use a Python binary with a version version >= 3.7  
```virtualenv -p=/usr/bin/python3.7 .my_venv```  
```source .my_venv/bin/activate```  
4. Update pip  
``` pip install --upgrade pip```  
5. Install dependencies  
```pip install -r requirements.txt```  

Tested on a clean install of ubuntu-20.10. System Requirements:
* git
* virtualenv

Starting with a clean installation of Ubuntu 20.10, code and
requirements will take up approximately 600MB.


## Running the simulation
The following requires the activation of the virtual environment you created during installation  
```source .my_venv/bin/activate```

I provide exemplary Jupyter Notebooks for [nursing homes](https://github.com/JanaLasser/agent_based_COVID_SEIRX/blob/dev/nursing_home/example_nursing_home.ipynb) and [schools](https://github.com/JanaLasser/agent_based_COVID_SEIRX/blob/dev/school/example_school.ipynb) that illustrate how a simulation model is set up and run for these two applications, how results are visualised and how data from a model run can be collected. Run the example notebook from the terminal:  
```jupyter-notebook example_nursing_home.ipynb```  

or  

```jupyter-notebook example_school.ipynb```

I also provide the [Jupyter Notebook](https://github.com/JanaLasser/agent_based_COVID_SEIRX/blob/dev/nursing_home/screening_frequency_data_creation.ipynb) used to run the simulations and create the data used in the publication **Agent-based simulations for optimized prevention of the spread of SARS-CoV-2 in nursing homes** as well as the [Jupyter Notebook](https://github.com/JanaLasser/agent_based_COVID_SEIRX/blob/master/nursing_home/screening_frequency_analysis.ipynb) used to create the heatmaps for the different analysed screnarios from the simulation data. Run these notebooks from the terminal using:

```jupyter-notebook nursing_home/screening_frequency_data_creation.ipynb```  

and  

```jupyter-notebook nursing_home/screening_frequency_analysis.ipynb```  


## Acknowledgements
I would like to thank Thomas Wochele-Thoma from [Caritas Austria](https://www.caritas.at/) for the fruitful discussions that led to the development of the nursing home scenario and for supplying the empirical outbreak observation data for this scenario. I would like to thank Lukas Richter and Daniela Schmid from [AGES Österreich](https://www.ages.at/startseite/) for supplying the data for the school scenario. I would like to than the numerous school directors and teachers that were willing to tell us about the daily life in Austrian schools during the pandemic. Lats but not least I would like to thank [Peter Klimek](https://www.csh.ac.at/researcher/peter-klimek/) from the [Complexity Science Hub Vienna](https://www.csh.ac.at/) and [Medical University Vienna](https://www.meduniwien.ac.at/web/) for the close collaboration during the inception and development of this project.

## References
Backer Jantien A, Klinkenberg Don, Wallinga Jacco. Incubation period of 2019 novel coronavirus (2019-nCoV) infections among travellers from Wuhan, China, 20–28 January 2020. Euro Surveill. 2020;25(5):pii=2000062. [DOI: 10.2807/1560-7917.ES.2020.25.5.2000062](https://doi.org/10.2807/1560-7917.ES.2020.25.5.2000062)

Ferretti, Luca, et al. "Quantifying SARS-CoV-2 transmission suggests epidemic control with digital contact tracing." Science 368.6491 (2020). [DOI: 10.1126/science.abb6936](https://doi.org/10.1126/science.abb6936)

He, X., Lau, E. H., Wu, P., Deng, X., Wang, J., Hao, X., ... & Mo, X. (2020). Temporal dynamics in viral shedding and transmissibility of COVID-19. Nature medicine, 26(5), 672-675. [DOI: 10.1038/s41591-020-0869-5](https://doi.org/10.1038/s41591-020-0869-5)  

Lauer, S. A., Grantz, K. H., Bi, Q., Jones, F. K., Zheng, Q., Meredith, H. R., ... & Lessler, J. (2020). The incubation period of coronavirus disease 2019 (COVID-19) from publicly reported confirmed cases: estimation and application. Annals of internal medicine, 172(9), 577-582. [DOI: 10.7326/M20-0504](https://doi.org/10.7326/M20-0504)  

Li, Q., Guan, X., Wu, P., Wang, X., Zhou, L., Tong, Y., ... & Xing, X. (2020). Early transmission dynamics in Wuhan, China, of novel coronavirus–infected pneumonia. New England Journal of Medicine. [DOI: 10.1056/NEJMoa2001316](https://doi.org/10.1056/NEJMoa2001316)  

Linton, N. M., Kobayashi, T., Yang, Y., Hayashi, K., Akhmetzhanov, A. R., Jung, S. M., ... & Nishiura, H. (2020). Incubation period and other epidemiological characteristics of 2019 novel coronavirus infections with right truncation: a statistical analysis of publicly available case data. Journal of clinical medicine, 9(2), 538. [DOI: 10.3390/jcm9020538](https://doi.org/10.3390/jcm9020538)  

Madewell, Z. J., Yang, Y., Longini, I. M., Halloran, M. E., & Dean, N. E. (2020). Household Transmission of SARS-CoV-2: A Systematic Review and Meta-analysis. JAMA network open, 3(12), [DOI10.1001/jamanetworkopen.2020.31756](https://doi.org/10.1001/jamanetworkopen.2020.31756).

McMichael, T. M., Currie, D. W., Clark, S., Pogosjans, S., Kay, M., Schwartz, N. G., ... & Ferro, J. (2020). Epidemiology of Covid-19 in a long-term care facility in King County, Washington. New England Journal of Medicine, 382(21), 2005-2011. [DOI: 10.1056/NEJMoa2005412](https://doi.org/10.1056/NEJMoa2005412) 

Nikolai, L. A., Meyer, C. G., Kremsner, P. G., & Velavan, T. P. (2020). Asymptomatic SARS Coronavirus 2 infection: Invisible yet invincible. International Journal of Infectious Diseases. [DOI: 10.1016/j.ijid.2020.08.076](https://doi.org/10.1016/j.ijid.2020.08.076)  

Pan, J., Harb, C., Leng, W., & Marr, L. C. (2020). Inward and outward effectiveness of cloth masks, a surgical mask, and a face shield. medRxiv. [DOI: 10.1101/2020.11.18.20233353](https://doi.org/10.1101/2020.11.18.20233353)

Polack, F.P., Thomas, S.J., Kitchin, N., Absalon, J., Gurtman, A., Lockhart, S., Perez, J.L., Pérez, M.G., et al. (2020). N Engl J Med. [DOI: 10.1056/NEJMoa2034577](10.1056/NEJMoa2034577)

Voysey, M., Clemens, S.A.C., Mahi, S.A., Weckx, L.Y., Folegatti, M.D., Aley, P.K., et al. (2020). Safety and efficacy of the ChAdOx1 nCoV-19 vaccine (AZD1222) against SARS-CoV-2: an interim analysis of four randomised controlled trials in Brazil, South Africa, and the UK. The Lancet. [DOI:10.1016/S0140-6736(20)32661-1](https://doi.org/10.1016/S0140-6736(20)32661-1)

Walsh, K. A., Jordan, K., Clyne, B., Rohde, D., Drummond, L., Byrne, P., ... & O'Neill, M. (2020). SARS-CoV-2 detection, viral load and infectivity over the course of an infection: SARS-CoV-2 detection, viral load and infectivity. Journal of Infection. [DOI: 10.1016/j.jinf.2020.06.067](10.1016/j.jinf.2020.06.067)  

Wu, J. T., Leung, K., & Leung, G. M. (2020). Nowcasting and forecasting the potential domestic and international spread of the 2019-nCoV outbreak originating in Wuhan, China: a modelling study. The Lancet, 395(10225), 689-697. [DOI: 10.1016/S0140-6736(20)30260-9](https://doi.org/10.1016/S0140-6736(20)30260-9) 

You, C., Deng, Y., Hu, W., Sun, J., Lin, Q., Zhou, F., ... & Zhou, X. H. (2020). Estimation of the time-varying reproduction number of COVID-19 outbreak in China. International Journal of Hygiene and Environmental Health, 113555. [DOI: 10.1016/j.ijheh.2020.113555](https://doi.org/10.1016/j.ijheh.2020.113555)
