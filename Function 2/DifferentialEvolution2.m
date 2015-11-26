function varargout = DifferentialEvolution2(varargin)
% DIFFERENTIALEVOLUTION2 MATLAB code for DifferentialEvolution2.fig
%      DIFFERENTIALEVOLUTION2, by itself, creates a new DIFFERENTIALEVOLUTION2 or raises the existing
%      singleton*.
%
%      H = DIFFERENTIALEVOLUTION2 returns the handle to a new DIFFERENTIALEVOLUTION2 or the handle to
%      the existing singleton*.
%
%      DIFFERENTIALEVOLUTION2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIFFERENTIALEVOLUTION2.M with the given input arguments.
%
%      DIFFERENTIALEVOLUTION2('Property','Value',...) creates a new DIFFERENTIALEVOLUTION2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DifferentialEvolution2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DifferentialEvolution2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DifferentialEvolution2

% Last Modified by GUIDE v2.5 26-Nov-2015 00:12:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DifferentialEvolution2_OpeningFcn, ...
                   'gui_OutputFcn',  @DifferentialEvolution2_OutputFcn, ...
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


% --- Executes just before DifferentialEvolution2 is made visible.
function DifferentialEvolution2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DifferentialEvolution2 (see VARARGIN)

% Choose default command line output for DifferentialEvolution2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DifferentialEvolution2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DifferentialEvolution2_OutputFcn(~, eventdata, handles) 
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
NumGenes = 10;
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
handles.axes1 = plot(XAxisData, YAxisData(1,:)); ylabel('Population Average Fitness'); xlabel('Iterations'); grid on; xlim([0,Generations]);
axes(handles.axes2);
handles.axes2 = plot(XAxisData, YAxisData(2,:)); ylabel('Best Individual Fitness'); xlabel('Iterations'); grid on; xlim([0,Generations]);


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
    [Individual, AverageFitness, MaxFit, MaxActualFit] = NewPopulation(Individual);    
end

dlmwrite('Function2DEResults.csv', [XAxisData; YAxisData],'-append', 'delimiter', '\t', 'newline', 'pc') 

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
   % Encode as two 5 bit genes. 1st & 6th bit determines sign. -15 -> 15
    TempX = bin2dec(char(Individual(k).Gene(2:5) + '0'));
    if Individual(k).Gene(1) == 1
        TempX = -TempX;        
    end
    
    TempY = bin2dec(char(Individual(k).Gene(7:10) + '0'));
    if Individual(k).Gene(6) == 1
        TempY = -TempY;        
    end
    
    ActualFit = 0.26*(TempX^2 + TempY^2) - 0.48*TempX*TempY;
    Individual(k).ActualFitness = ActualFit;
    
    % Convert to a maximisation problem
    TempFit = (1/(1+ActualFit));
 
    % Write fitness to individual
    Individual(k).Fitness = TempFit;
    
    % For use in working out individual probability
    TotalFitness = TotalFitness + TempFit;
end


function [ Individual, AverageFitness, MaxFit, MaxActualFit ] = NewPopulation(Individual)
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

% Run through "population" number of times before plotting
for k = 1:1
    
    % Calculate individual probability
    for m = 1:Population
        Individual(m).Probability = Individual(m).Fitness / TotalFitness;
    end

    % Create "roulette wheel"
    SelectionWheel = cumsum([Individual.Probability]);
      
    % 3rd solution picked based on fitness
    Selection = rand;
    for n = 1:Population
       if SelectionWheel(n) >= Selection
          MutationInput(1:4) = n; 
          break;
       end
    end
    
    % Differentials picked at random
    while (MutationInput(2)== MutationInput(1))
        MutationInput(2) = round(rand*(Population-1))+1;
    end
    
    while (MutationInput(3)== MutationInput(2)) || (MutationInput(3)== MutationInput(1))
        MutationInput(3) = round(rand*(Population-1))+1;
    end
    
    % decode x + y for all 3 individuals
    for q = 1:3
        DecodeVal(q,1) = bin2dec(char(Individual(MutationInput(q)).Gene(2:5) + '0'));
        if Individual(MutationInput(q)).Gene(1) == 1
            DecodeVal(q,1) = -DecodeVal(q,1);        
        end

        DecodeVal(q,2) = bin2dec(char(Individual(MutationInput(q)).Gene(7:10) + '0'));
        if Individual(MutationInput(q)).Gene(6) == 1
            DecodeVal(q,2) = -DecodeVal(q,2);        
        end
    end
    
    % add to 3rd variable
    DonorVal(1) = round(DecodeVal(1,1) + MutationRate*(DecodeVal(2,1) - DecodeVal(3,1)));
    DonorVal(2) = round(DecodeVal(1,2) + MutationRate*(DecodeVal(2,2) - DecodeVal(3,2)));
    
    % Make sure new values don't exceed limits
    DonorVal(DonorVal>15) = 15;
    DonorVal(DonorVal<-15) = -15;

    % Rencode
    for q = 1:2
        if DonorVal(q) < 0
            Donor.Gene(-3+(5*q):(5*q)) = dec2bin(-DonorVal(q),4)=='1';
            Donor.Gene(-4+(5*q)) = 1;
        else
            Donor.Gene(-3+(5*q):(5*q)) = dec2bin(DonorVal(q),4)=='1';
            Donor.Gene(-4+(5*q)) = 0;
        end
    end
    
    % Target solution picked based on fitness
    while (MutationInput(4)== MutationInput(3)) || (MutationInput(4)== MutationInput(2)) || (MutationInput(4)== MutationInput(1))
        Selection = rand;
        for n = 1:Population
           if SelectionWheel(n) >= Selection
              MutationInput(4) = n; 
              break;
           end
        end
    end
    
    % Perform Recombination - Binomial crossover
    Irand = round(rand*(NumGenes-1))+1;
    for q = 1:NumGenes
        if (rand > CrossOverRate) && (q ~= Irand)
           Donor.Gene(q) =  Individual(MutationInput(4)).Gene(q);
        end   
    end
    
    
    Donor.ActualFitness = 0; 
    
    % Calculate fitness using fitness function
    DonorVal(1) = bin2dec(char(Donor.Gene(2:5) + '0'));
    if Donor.Gene(1) == 1
        DonorVal(1) = -DonorVal(1);        
    end
    DonorVal(2) = bin2dec(char(Donor.Gene(7:10) + '0'));
    if Donor.Gene(6) == 1
        DonorVal(2) = -DonorVal(2);        
    end
    
    ActualFit = 0.26*(DonorVal(1)^2 + DonorVal(2)^2) - 0.48*DonorVal(1)*DonorVal(2);
    Donor.ActualFitness = ActualFit;
    
    % Convert to a maximisation problem
    TempFit = (1/(1+ActualFit));
    % Write fitness to individual
    Donor.Fitness = TempFit;
    
    % Compare Donor and Target individual, keep the best one
    if Donor.Fitness > Individual(MutationInput(4)).Fitness
        Individual(MutationInput(4)).Gene = Donor.Gene;
    end


    % Calculate Offspring Fitness
    [Individual] = CalculateFitness(Individual);

end

% Work out Maxfit & average fitness for plotting
[MaxFit, p] = max([Individual.Fitness]);
MaxActualFit = Individual(p).ActualFitness;
AverageFitness = mean([Individual.ActualFitness]);



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
