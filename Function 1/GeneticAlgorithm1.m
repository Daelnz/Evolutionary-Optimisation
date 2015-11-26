function varargout = GeneticAlgorithm1(varargin)
% GENETICALGORITHM1 MATLAB code for GeneticAlgorithm1.fig
%      GENETICALGORITHM1, by itself, creates a new GENETICALGORITHM1 or raises the existing
%      singleton*.
%
%      H = GENETICALGORITHM1 returns the handle to a new GENETICALGORITHM1 or the handle to
%      the existing singleton*.
%
%      GENETICALGORITHM1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GENETICALGORITHM1.M with the given input arguments.
%
%      GENETICALGORITHM1('Property','Value',...) creates a new GENETICALGORITHM1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GeneticAlgorithm1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GeneticAlgorithm1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GeneticAlgorithm1

% Last Modified by GUIDE v2.5 24-Nov-2015 21:55:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GeneticAlgorithm1_OpeningFcn, ...
                   'gui_OutputFcn',  @GeneticAlgorithm1_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GeneticAlgorithm1 is made visible.
function GeneticAlgorithm1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GeneticAlgorithm1 (see VARARGIN)

% Choose default command line output for GeneticAlgorithm1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GeneticAlgorithm1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GeneticAlgorithm1_OutputFcn(~, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in StartGA.
function StartGA_Callback(hObject, eventdata, handles)
% hObject    handle to StartGA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Population;
global NumGenes;
global MutationRate;
global CrossOverRate;
global Generations;

% Get values from GUI
Population = str2double(get(handles.PopulationSize, 'String'));
NumGenes = 8;
MutationRate = str2double(get(handles.MutationProb, 'String'));
CrossOverRate = str2double(get(handles.CrossoverProb, 'String'));
Generations = str2double(get(handles.NumGenerations, 'String'));

clc;

% Set graph to be length of amount of generations
XAxisData = linspace(1,Generations,Generations);

% Populate Data with zeroes to start with
YAxisData = zeros(2,Generations);

% Setups graph details
axes(handles.axes1); %Ensures correct graph is selected
handles.axes1 = plot(XAxisData, YAxisData(1,:)); ylabel('Population Average Fitness'); xlabel('Generations'); grid on; xlim([0,Generations]);
axes(handles.axes2);
handles.axes2 = plot(XAxisData, YAxisData(2,:)); ylabel('Best Individual Fitness'); xlabel('Generations'); grid on; xlim([0,Generations]);


%% Genetic Algorithm Code

% Generate empty structure for population genes and fitness
Field1 = 'Gene';
Field2 = 'Fitness';
Field3 = 'Probability';
Field4 = 'ActualFitness';

Individual = struct(Field1, [], Field2, [], Field3, [], Field4, []);

% Generate population with random fitness
for j = 1:Population
    for i = 1:NumGenes
        Individual(j).Gene(i) = round(rand);
    end
end

% Calculate Fitness
[Individual] = CalculateFitness(Individual);

% Work out Maxfit & average fitness for plotting
[MaxFit, p] = max([Individual.Fitness]);
MaxActualFit = Individual(p).ActualFitness;
AverageFitness = mean([Individual.ActualFitness]);

% Loop through egenerations
for t = 1:Generations
    
    % Plot Total Fitness + Max fitness
    YAxisData(1,t) = AverageFitness;
    YAxisData(2,t) = MaxActualFit;  
    set(handles.axes1, 'xdata', XAxisData, 'ydata', YAxisData(1,:));
    set(handles.axes2, 'xdata', XAxisData, 'ydata', YAxisData(2,:));
    pause(0.01); % Allow plot to update
    
    % Generate Offspring + Calculate Fitness
    [Individual, AverageFitness, MaxFit, MaxActualFit] = NewPopulation(Individual, MaxFit);    
end

dlmwrite('Function1Results.csv', [XAxisData; YAxisData],'-append', 'delimiter', '\t', 'newline', 'pc')  

% Work out how many generations it took to find the best Individual
for b = 1:Generations
    if (YAxisData(2,b) == MaxActualFit)
        set(handles.BestGeneration, 'String', num2str(b));
       break;
    end
end

% Find best actual fitness
for k = 1:Population
    if (Individual(k).Fitness == MaxFit)
        set(handles.BestFit, 'String', Individual(k).ActualFitness);
        break;
    end
end



function [ Individual] = CalculateFitness(Individual)
%CalculateFitness Calculates the fitness of each individual
%   Counts the number of 1's in each indiividuals genes to calculate it's
%   fitness. Also returns AverageFitness of the population

global Population;
global NumGenes;
global TotalFitness;

TotalFitness = 0;

% Calculate Fitness of each individual
for k = 1:Population
    
    Individual(k).ActualFitness = 0; 
    
    % Calculate fitness using fitness function
    TempFit = bin2dec(char(Individual(k).Gene + '0')); 
    TempFit = (TempFit^2);
 
    % Write fitness to individual
    Individual(k).Fitness = TempFit;
    Individual(k).ActualFitness = TempFit;
    
    % For use in working out individual probability
    TotalFitness = TotalFitness + TempFit; 
end


function [ Offspring, AverageFitness, MaxFit, MaxActualFit ] = NewPopulation(Individual, MaxFit)
%NewPopulation generates a new population of individuals from the parent
%population
%   The function randomly selects parents, uses their fitness to determine
%   the childs genes then randomly applies crossover or mutation based on
%   the probability of each.

global Population;
global NumGenes;
global MutationRate;
global CrossOverRate;
global TotalFitness;

% Calaculate individual probability
for m = 1:Population
    Individual(m).Probability = Individual(m).Fitness / TotalFitness;
end

% Create "roulette wheel"
SelectionWheel = cumsum([Individual.Probability]);

% Find best individual from parent population
[~, q] = max([Individual(:).Fitness]);
BestIndividual = Individual(q);

% Create New Population
for k = 1:Population
      
    % Pick Parent 1
    Selection = rand;
    for n = 1:Population
       if SelectionWheel(n) >= Selection
          Parent1 = n; 
          break;
       end
    end
    
    % Pick Parent 2
    Selection = rand;
    for n = 1:Population
       if SelectionWheel(n) >= Selection
          Parent2 = n; 
          break;
       end
    end

        % Apply Crossover
        CrossOverPoint = round(rand*(NumGenes-1))+1; 
        
        % Offspring based on highest fitness parent + probability of crossover 
        if Individual(Parent1).Fitness >= Individual(Parent2).Fitness
            Offspring(k).Gene = Individual(Parent1).Gene;
            if rand < CrossOverRate
                Offspring(k).Gene(CrossOverPoint:end) = Individual(Parent2).Gene(CrossOverPoint:end);
            end
        else
            Offspring(k).Gene = Individual(Parent2).Gene;
            if rand < CrossOverRate
                Offspring(k).Gene(CrossOverPoint:end) = Individual(Parent1).Gene(CrossOverPoint:end);
            end
        end
        
        % Impliment Mutation
        for l = 1:NumGenes
             if rand < MutationRate
                Offspring(k).Gene(l) = ~Offspring(k).Gene(l);
             end
        end
end

% Calculate Offspring Fitness
[Offspring] = CalculateFitness(Offspring);

% Find worst new individual and replace with best individual from parents generation
[~, MinIndividual] = min([Offspring.Fitness]);
TotalFitness = TotalFitness - Offspring(MinIndividual).Fitness + BestIndividual.Fitness; 
Offspring(MinIndividual).Gene = BestIndividual.Gene;
Offspring(MinIndividual).Fitness = BestIndividual.Fitness;
Offspring(MinIndividual).ActualFitness = BestIndividual.ActualFitness;


% Work out Maxfit & average fitness for plotting
[MaxFit, p] = max([Offspring.Fitness]);
MaxActualFit = Offspring(p).ActualFitness;
AverageFitness = mean([Offspring.ActualFitness]);



function PopulationSize_Callback(hObject, eventdata, handles)
% hObject    handle to PopulationSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PopulationSize as text
%        str2double(get(hObject,'String')) returns contents of PopulationSize as a double


% --- Executes during object creation, after setting all properties.
function PopulationSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PopulationSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NSize_Callback(hObject, eventdata, handles)
% hObject    handle to NSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NSize as text
%        str2double(get(hObject,'String')) returns contents of NSize as a double


% --- Executes during object creation, after setting all properties.
function NSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NumGenerations_Callback(hObject, eventdata, handles)
% hObject    handle to NumGenerations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumGenerations as text
%        str2double(get(hObject,'String')) returns contents of NumGenerations as a double


% --- Executes during object creation, after setting all properties.
function NumGenerations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumGenerations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MutationProb_Callback(hObject, eventdata, handles)
% hObject    handle to MutationProb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MutationProb as text
%        str2double(get(hObject,'String')) returns contents of MutationProb as a double


% --- Executes during object creation, after setting all properties.
function MutationProb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MutationProb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CrossoverProb_Callback(hObject, eventdata, handles)
% hObject    handle to CrossoverProb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CrossoverProb as text
%        str2double(get(hObject,'String')) returns contents of CrossoverProb as a double


% --- Executes during object creation, after setting all properties.
function CrossoverProb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CrossoverProb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BestGeneration_Callback(hObject, eventdata, handles)
% hObject    handle to BestGeneration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BestGeneration as text
%        str2double(get(hObject,'String')) returns contents of BestGeneration as a double


% --- Executes during object creation, after setting all properties.
function BestGeneration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BestGeneration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BestFit_Callback(hObject, eventdata, handles)
% hObject    handle to BestFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BestFit as text
%        str2double(get(hObject,'String')) returns contents of BestFit as a double


% --- Executes during object creation, after setting all properties.
function BestFit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BestFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
