import numpy as np
import networkx as nx
from mesa import Model
from mesa.time import RandomActivation, SimultaneousActivation
from mesa.datacollection import DataCollector

import sys
sys.path.insert(0,'..')
from model_SEIRX import*


## data collection functions ##

def count_E_student(model):
    E = np.asarray(
        [a.exposed for a in model.schedule.agents if a.type == 'student']).sum()
    return E


def count_I_student(model):
    I = np.asarray(
        [a.infectious for a in model.schedule.agents if a.type == 'student']).sum()
    return I


def count_I_symptomatic_student(model):
    I = np.asarray([a.infectious for a in model.schedule.agents if
        (a.type == 'student'and a.symptomatic_course)]).sum()
    return I


def count_I_asymptomatic_student(model):
    I = np.asarray([a.infectious for a in model.schedule.agents if
        (a.type == 'student'and a.symptomatic_course == False)]).sum()
    return I


def count_R_student(model):
    R = np.asarray(
        [a.recovered for a in model.schedule.agents if a.type == 'student']).sum()
    return R


def count_X_student(model):
    X = np.asarray(
        [a.quarantined for a in model.schedule.agents if a.type == 'student']).sum()
    return X


def count_E_teacher(model):
    E = np.asarray(
        [a.exposed for a in model.schedule.agents if a.type == 'teacher']).sum()
    return E


def count_I_teacher(model):
    I = np.asarray(
        [a.infectious for a in model.schedule.agents if a.type == 'teacher']).sum()
    return I


def count_I_symptomatic_teacher(model):
    I = np.asarray([a.infectious for a in model.schedule.agents if
        (a.type == 'teacher'and a.symptomatic_course)]).sum()
    return I


def count_I_asymptomatic_teacher(model):
    I = np.asarray([a.infectious for a in model.schedule.agents if
        (a.type == 'teacher'and a.symptomatic_course == False)]).sum()
    return I


def count_R_teacher(model):
    R = np.asarray(
        [a.recovered for a in model.schedule.agents if a.type == 'teacher']).sum()
    return R


def count_X_teacher(model):
    X = np.asarray(
        [a.quarantined for a in model.schedule.agents if a.type == 'teacher']).sum()
    return X


def count_E_family_member(model):
    E = np.asarray(
        [a.exposed for a in model.schedule.agents if a.type == 'family_member']).sum()
    return E


def count_I_family_member(model):
    I = np.asarray(
        [a.infectious for a in model.schedule.agents if a.type == 'family_member']).sum()
    return I


def count_I_symptomatic_family_member(model):
    I = np.asarray([a.infectious for a in model.schedule.agents if
        (a.type == 'family_member'and a.symptomatic_course)]).sum()
    return I


def count_I_asymptomatic_family_member(model):
    I = np.asarray([a.infectious for a in model.schedule.agents if
        (a.type == 'family_member'and a.symptomatic_course == False)]).sum()
    return I


def count_R_family_member(model):
    R = np.asarray(
        [a.recovered for a in model.schedule.agents if a.type == 'family_member']).sum()
    return R


def count_X_family_member(model):
    X = np.asarray(
        [a.quarantined for a in model.schedule.agents if a.type == 'family_member']).sum()
    return X


def check_reactive_student_screen(model):
    return model.screened_agents['reactive']['student']


def check_follow_up_student_screen(model):
    return model.screened_agents['follow_up']['student']


def check_preventive_student_screen(model):
    return model.screened_agents['preventive']['student']


def check_reactive_teacher_screen(model):
    return model.screened_agents['reactive']['teacher']


def check_follow_up_teacher_screen(model):
    return model.screened_agents['follow_up']['teacher']


def check_preventive_teacher_screen(model):
    return model.screened_agents['preventive']['teacher']


def check_reactive_family_member_screen(model):
    return model.screened_agents['reactive']['family_member']


def check_follow_up_family_member_screen(model):
    return model.screened_agents['follow_up']['family_member']


def check_preventive_family_member_screen(model):
    return model.screened_agents['preventive']['family_member']



data_collection_functions = \
    {
    'student':
        {
        'E':count_E_student,
        'I':count_I_student,
        'I_asymptomatic':count_I_asymptomatic_student,
        'I_symptomatic':count_I_symptomatic_student,
        'R':count_R_student,
        'X':count_X_student
         },
    'teacher':
        {
        'E':count_E_teacher,
        'I':count_I_teacher,
        'I_asymptomatic':count_I_asymptomatic_teacher,
        'I_symptomatic':count_I_symptomatic_teacher,
        'R':count_R_teacher,
        'X':count_X_teacher
         },
    'family_member':
        {
        'E':count_E_family_member,
        'I':count_I_family_member,
        'I_asymptomatic':count_I_asymptomatic_family_member,
        'I_symptomatic':count_I_symptomatic_family_member,
        'R':count_R_family_member,
        'X':count_X_family_member
         }
    }



class SEIRX_school(SEIRX):
    '''
    Model specific parameters:
        age_risk_discount: discount factor that lowers the transmission and 
        reception risk of agents based on age for children. This is only applied
        to student agents as all other agents are assumed to be adults. This 
        parameter needs to be calibrated against data.

    See documentation of model_SEIRX for the description of other parameters.
    '''

    def __init__(self, G, verbosity=0, 
        base_transmission_risk = 0.05, testing=True,
        exposure_duration=4, time_until_symptoms=6, infection_duration=11, 
        quarantine_duration=14, subclinical_modifier=1,
        infection_risk_contact_type_weights={
            'very_far': 0.1, 'far': 0.25, 'intermediate': 0.5, 'close': 1},
        K1_contact_types=['close'], diagnostic_test_type='one_day_PCR',
        preventive_screening_test_type='one_day_PCR',
        follow_up_testing_interval=None, liberating_testing=False,
        index_case='employee', 
        agent_types={'type1': {'screening_interval': None,
                              'index_probability': None,
                              'mask':False}},
        age_transmission_risk_discount = {'slope':-0.05, 'intercept':1},
        age_symptom_discount = {'slope':-0.02545, 'intercept':0.854545},
        mask_filter_efficiency = {'exhale':0, 'inhale':0},
        transmission_risk_ventilation_modifier = 0,
        seed=None):


        super().__init__(G, verbosity, base_transmission_risk, testing,
            exposure_duration, time_until_symptoms, infection_duration,
            quarantine_duration, subclinical_modifier,
            infection_risk_contact_type_weights,
            K1_contact_types, diagnostic_test_type,
            preventive_screening_test_type,
            follow_up_testing_interval, liberating_testing,
            index_case, agent_types, age_transmission_risk_discount,
            age_symptom_discount, mask_filter_efficiency, 
            transmission_risk_ventilation_modifier, seed)
        
        self.MG = G
        self.weekday_connections = {}
        all_edges = self.MG.edges(keys=True, data='weekday')
        N_weekdays = 7
        for i in range(1, N_weekdays + 1):
            wd_edges = [(u, v, k) for (u, v, k, wd) in all_edges if wd == i]
            self.weekday_connections[i] = G.edge_subgraph(wd_edges).copy()
        
        
        # data collectors to save population counts and agent states every
        # time step
        model_reporters = {}
        for agent_type in self.agent_types:
            for state in ['E','I','I_asymptomatic','I_symptomatic','R','X']:
                model_reporters.update({'{}_{}'.format(state, agent_type):\
                    data_collection_functions[agent_type][state]})

        model_reporters.update(\
            {
            'screen_students_reactive':check_reactive_student_screen,
            'screen_students_follow_up':check_follow_up_student_screen,
            'screen_students_preventive':check_preventive_student_screen,
            'screen_teachers_reactive':check_reactive_teacher_screen,
            'screen_teachers_follow_up':check_follow_up_teacher_screen,
            'screen_teachers_preventive':check_preventive_teacher_screen,
            'screen_family_members_reactive':check_reactive_family_member_screen,
            'screen_family_members_follow_up':check_follow_up_family_member_screen,
            'screen_family_members_preventive':check_preventive_family_member_screen,
            'N_diagnostic_tests':get_N_diagnostic_tests,
            'N_preventive_screening_tests':get_N_preventive_screening_tests,
            'undetected_infections':get_undetected_infections,
            'predetected_infections':get_predetected_infections,
            'pending_test_infections':get_pending_test_infections
            })

        agent_reporters =\
            {
            'infection_state':get_infection_state,
            'quarantine_state':get_quarantine_state
             }

        self.datacollector = DataCollector(
            model_reporters = model_reporters,
            agent_reporters = agent_reporters)

    def calculate_transmission_probability(self, source, target, base_risk):
        """
        Calculates the risk of transmitting an infection between a source agent
        and a target agent given the model's and agent's properties and the base
        transmission risk.

        Transmission is an independent Bernoulli trial with a probability of
        success p. The probability of transmission without any modifications
        by for example masks or ventilation is given by the base_risk, which
        is calibrated in the model. The probability is modified by contact type
        q1 (also calibrated in the model), age of the transmitting agent q2 
        & age of the receiving agent q3 (both age dependencies are linear in 
        age and the same, and they are calibrated), infection progression q4
        (from literature), reduction of exhaled viral load of the source by mask
        wearing q5 (from literature), reduction of inhaled viral load by the
        target q6 (from literature), and ventilation of the rooms q7 (from
        literature).

        Parameters
        ----------
        source : agent_SEIRX
            Source agent that transmits the infection to the target.
        target: agent_SEIRX
            Target agent that (potentially) receives the infection from the 
            source.
        base_risk : float
            Probability p of infection transmission without any modifications
            through prevention measures.

        Returns
        -------
        p : float
            Modified transmission risk.
        """
        n1 = source.ID
        n2 = target.ID
        tmp = [n1, n2]
        tmp.sort()
        n1, n2 = tmp
        key = n1 + n2 + 'd{}'.format(self.weekday)
        link_type = self.G.get_edge_data(n1, n2, key)['link_type']

        q1 = self.get_transmission_risk_contact_type_modifier(source, target)
        q2 = self.get_transmission_risk_age_modifier_transmission(source)
        q3 = self.get_transmission_risk_age_modifier_reception(target)
        q4 = self.get_transmission_risk_progression_modifier(source)
        q5 = self.get_transmission_risk_subclinical_modifier(source)

        # contact types where masks and ventilation are irrelevant
        if link_type in ['student_household', 'teacher_household']:
            p = 1 - (1 - base_risk * (1- q1) * (1 - q2) * (1 - q3) * \
                (1 - q4) * (1 - q5))

        # contact types were masks and ventilation are relevant
        elif link_type in ['student_student_intra_class',
                           'student_student_table_neighbour',
                           'student_student_daycare',
                           'teacher_teacher_short',
                           'teacher_teacher_long',
                           'teacher_teacher_team_teaching',
                           'teacher_teacher_daycare_supervision',
                           'teaching_teacher_student',
                           'daycare_supervision_teacher_student']:
            q6 = self.get_transmission_risk_exhale_modifier(source)
            q7 = self.get_transmission_risk_inhale_modifier(target)
            q8 = self.get_transmission_risk_ventilation_modifier()

            p = 1 - (1 - base_risk * (1- q1) * (1 - q2) * (1 - q3) * \
                (1 - q4) * (1 - q5) * (1 - q6) * (1 - q7) * (1 - q8))

        else:
            print('unknown link type: {}'.format(link_type))
            p = None
        return p


    def step(self):
        self.weekday = (self.Nstep + self.weekday_offset)% 7 + 1
        self.G = self.weekday_connections[self.weekday]
        if self.verbosity > 0:
            print('weekday {}'.format(self.weekday))

        if self.testing:
            for agent_type in self.agent_types:
                for screen_type in ['reactive', 'follow_up', 'preventive']:
                    self.screened_agents[screen_type][agent_type] = False

            if self.verbosity > 0: 
                print('* testing and tracing *')
            
            self.test_symptomatic_agents()
            

            # collect and act on new test results
            agents_with_test_results = self.collect_test_results()
            for a in agents_with_test_results:
                a.act_on_test_result()
            
            self.quarantine_contacts()

            # screening:
            # a screen should take place if
            # (a) there are new positive test results
            # (b) as a follow-up screen for a screen that was initiated because
            # of new positive cases
            # (c) if there is a preventive screening policy and it is time for
            # a preventive screen in a given agent group

            # (a)
            if (self.testing == 'background' or self.testing == 'preventive')\
               and self.new_positive_tests == True:
                for agent_type in ['teacher', 'student']:
	                self.screen_agents(
	                    agent_type, self.Testing.diagnostic_test_type, 'reactive')
	                self.scheduled_follow_up_screen[agent_type] = True

            # (b)
            elif (self.testing == 'background' or self.testing == 'preventive') and \
                self.Testing.follow_up_testing_interval != None and \
                sum(list(self.scheduled_follow_up_screen.values())) > 0:
                for agent_type in ['teacher', 'student']:
                    if self.scheduled_follow_up_screen[agent_type] and\
                       self.days_since_last_agent_screen[agent_type] >=\
                       self.Testing.follow_up_testing_interval:
                        self.screen_agents(
                            agent_type, self.Testing.diagnostic_test_type, 'follow_up')
                    else:
                        if self.verbosity > 0: 
                            print('not initiating {} follow-up screen (last screen too close)'\
                            	.format(agent_type))

            # (c) 
            elif self.testing == 'preventive' and \
                np.any(list(self.Testing.screening_intervals.values())):

                for agent_type in ['teacher', 'student']:
                    interval = self.Testing.screening_intervals[agent_type]
                    assert interval in [7, 3, 2, None], \
                        'testing interval {} for agent type {} not supported!'\
                        .format(interval, agent_type)

                    # (c.1) testing every 7 days = testing on mondays
                    if interval == 7 and self.weekday == 1:
                        self.screen_agents(agent_type,
                            self.Testing.preventive_screening_test_type, 'preventive')
                    # (c.2) testing every 3 days = testing on mondays & thursdays
                    elif interval == 3 and self.weekday in [1, 4]:
                            self.screen_agents(agent_type,
                            self.Testing.preventive_screening_test_type, 'preventive')
                    # (c.3) testing every 2 days = testing on mondays, wednesdays & fridays
                    elif interval == 2 and self.weekday in [1, 3, 5]:
                            self.screen_agents(agent_type,
                            self.Testing.preventive_screening_test_type, 'preventive')
                    # No interval specified = no testing, even if testing mode == preventive
                    elif interval == None:
                        pass
                    else:
                        if self.verbosity > 0:
                            print('not initiating {} preventive screen (wrong weekday)'\
                                    .format(agent_type))
            else:
                # do nothing
                pass

            for agent_type in self.agent_types:
            	if not (self.screened_agents['reactive'][agent_type] or \
            		    self.screened_agents['follow_up'][agent_type] or \
            		    self.screened_agents['preventive'][agent_type]):
            			self.days_since_last_agent_screen[agent_type] += 1


        if self.verbosity > 0: print('* agent interaction *')
        self.datacollector.collect(self)
        self.schedule.step()
        self.Nstep += 1
