%[text] # Build a cranked wing from its CST parameters
%[text] Move a slider and the wing is rebuilt with `CST2Wing` (the second section runs again).
%[text] Kulfan (2008, Sec. XII-XIII) "scalar loft": every CST coefficient of the sections and every planform quantity is a Bernstein polynomial along the span `eta = y/b`. Here the span has one break (a kink at `eta_break`) and the polynomials are linear (`ny = ny_plan = 1`) in each of the two segments, so every control value is simply the value at the root, at the break or at the tip:
%[text] Default values: a transport-like wing, NACA 2415 at the root and NACA 2410 at the tip.
%[text] ## Setup
clc; clear; close all
addpath(fullfile('..','src'))       % matlab_cst functions (run from the examples folder)

% plot settings
axFS = 13;
axFN = "Times New Roman";
lw   = 2;
axlw = 1.25;
set(0, 'defaultAxesTickLabelInterpreter','latex', ...
       'defaultTextInterpreter','latex', ...
       'DefaultLegendInterpreter','latex');
%%
%[text] ## CST parameters
% span and break
b         = 18; %[control:slider:12d2]{"position":[13,15]}
eta_break = 0.35; %[control:slider:4ca8]{"position":[13,17]}
% chord at root, break, tip [m]
c_root  = 7; %[control:slider:bc99]{"position":[11,12]}
c_break = 4.5; %[control:slider:7004]{"position":[11,14]}
c_tip   = 1.6; %[control:slider:111b]{"position":[11,14]}
% untwisted leading edge at break and tip [m] (root at 0): sweep and dihedral
Xle_break = 3; %[control:slider:5587]{"position":[13,14]}
Xle_tip   = 9.5; %[control:slider:ef24]{"position":[13,16]}
Zle_break = 0.3; %[control:slider:3be2]{"position":[13,16]}
Zle_tip   = 1.2; %[control:slider:bb63]{"position":[13,16]}
% twist at root, break, tip [deg]
twist_root  = 3; %[control:slider:4f15]{"position":[15,16]}
twist_break = 1; %[control:slider:32ef]{"position":[15,16]}
twist_tip   = -1.5; %[control:slider:bc80]{"position":[15,19]}
% root airfoil (Al(1) = -Au(1))
AuR = zeros(4,1); AlR = zeros(4,1);
AuR(1) = 0.213; %[control:slider:e35f]{"position":[10,15]}
AuR(2) = 0.281; %[control:slider:ac7f]{"position":[10,15]}
AuR(3) = 0.197; %[control:slider:43f6]{"position":[10,15]}
AuR(4) = 0.268; %[control:slider:5127]{"position":[10,15]}
AlR(2) = -0.101; %[control:slider:8529]{"position":[10,16]}
AlR(3) = -0.156; %[control:slider:8ce2]{"position":[10,16]}
AlR(4) = -0.09; %[control:slider:c1ea]{"position":[10,15]}
% tip airfoil (Al(1) = -Au(1))
AuT = zeros(4,1); AlT = zeros(4,1);
AuT(1) = 0.142; %[control:slider:ecd7]{"position":[10,15]}
AuT(2) = 0.225; %[control:slider:0ce0]{"position":[10,15]}
AuT(3) = 0.134; %[control:slider:f37a]{"position":[10,15]}
AuT(4) = 0.213; %[control:slider:0b54]{"position":[10,15]}
AlT(2) = -0.029; %[control:slider:83d6]{"position":[10,16]}
AlT(3) = -0.101; %[control:slider:2f5c]{"position":[10,16]}
AlT(4) = -0.027; %[control:slider:1c26]{"position":[10,16]}

% BUILD ---------------------------------------------------------------------
AlR(1) = -AuR(1);  AlT(1) = -AuT(1);                     % same LE radius on both sides
AuB    = AuR + eta_break*(AuT - AuR);                     % break airfoil: linear root -> tip
AlB    = AlR + eta_break*(AlT - AlR);
CST            = struct();
CST.b          = b;                                       % m
CST.eta_breaks = eta_break;
CST.Bu         = [AuR'; AuB'; AuT'];                      % rows: root, break, tip
CST.Bl         = [AlR'; AlB'; AlT'];
CST.c          = [c_root; c_break; c_tip];                % m
CST.Xle        = [0; Xle_break; Xle_tip];                 % m
CST.Zle        = [0; Zle_break; Zle_tip];                 % m
CST.twist      = deg2rad([twist_root; twist_break; twist_tip]);   % rad
W = CST2Wing(CST,eta=unique([linspace(0,1,25)'; eta_break]),npsi=61);
fprintf('S = %.2f m^2 (half wing) | AR = %.2f | mac = %.2f m | LE sweep %.1f / %.1f deg | dihedral %.1f / %.1f deg\n', ... %[output:group:0fec2898] %[output:94a6dad1]
        W.S,W.AR,W.mac,rad2deg(W.sweepLE([1 end])),rad2deg(W.dihedral([1 end]))); %[output:group:0fec2898] %[output:94a6dad1]

% PLOT ----------------------------------------------------------------------
figure('Position',[100 100 1000 800]); %[output:818cf68e]
tiledlayout(2,2,'TileSpacing','compact'); %[output:818cf68e]
ax1 = nexttile([1 2]); %[output:818cf68e]
W.plot(ax1,mirror=true); %[output:818cf68e]
view(ax1,[-140 25]); axis(ax1,'off'); %[output:818cf68e]
ax2 = nexttile; hold on %[output:818cf68e]
W.plotPlanform(ax2,'k-','LineWidth',lw); %[output:818cf68e]
xline(ax2,W.b*eta_break,'r--','break'); %[output:818cf68e]
styleAx(ax2,axFS,axFN,axlw); %[output:818cf68e]
ax3 = nexttile; hold on %[output:818cf68e]
iSec = [1 find(abs(W.eta - eta_break) < 1e-9) numel(W.eta)];
names = ["root","break","tip"];
for k = 1:3
    W.Airfoil{iSec(k)}.plot(ax3,'-','LineWidth',lw,'DisplayName',names(k)); %[output:818cf68e]
end
xlabel('$x/c$'); ylabel('$z/c$'); legend('Location','northeast'); %[output:818cf68e]
styleAx(ax3,axFS,axFN,axlw); %[output:818cf68e]
%%
%[text] ## Local functions
function styleAx(ax,axFS,axFN,axlw)
% axes style of the script skeleton
ax.FontName  = axFN;
ax.FontSize  = axFS;
ax.LineWidth = axlw;
ax.XMinorTick = "on"; ax.YMinorTick = "on";
ax.Layer = "top";
box(ax,'on'); grid(ax,'on');
ax.GridAlpha = 0.1;
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
%[control:slider:12d2]
%   data: {"defaultValue":18,"label":"b","max":40,"min":5,"run":"Section","runOn":"ValueChanged","step":0.5}
%---
%[control:slider:4ca8]
%   data: {"defaultValue":0.35,"label":"eta_break","max":0.9,"min":0.1,"run":"Section","runOn":"ValueChanged","step":0.01}
%---
%[control:slider:bc99]
%   data: {"defaultValue":7,"label":"c_root","max":12,"min":1,"run":"Section","runOn":"ValueChanged","step":0.1}
%---
%[control:slider:7004]
%   data: {"defaultValue":4.5,"label":"c_break","max":10,"min":0.5,"run":"Section","runOn":"ValueChanged","step":0.1}
%---
%[control:slider:111b]
%   data: {"defaultValue":1.6,"label":"c_tip","max":6,"min":0.3,"run":"Section","runOn":"ValueChanged","step":0.1}
%---
%[control:slider:5587]
%   data: {"defaultValue":3,"label":"Xle_break","max":15,"min":-5,"run":"Section","runOn":"ValueChanged","step":0.1}
%---
%[control:slider:ef24]
%   data: {"defaultValue":9.5,"label":"Xle_tip","max":30,"min":-5,"run":"Section","runOn":"ValueChanged","step":0.1}
%---
%[control:slider:3be2]
%   data: {"defaultValue":0.3,"label":"Zle_break","max":3,"min":-3,"run":"Section","runOn":"ValueChanged","step":0.05}
%---
%[control:slider:bb63]
%   data: {"defaultValue":1.2,"label":"Zle_tip","max":8,"min":-5,"run":"Section","runOn":"ValueChanged","step":0.05}
%---
%[control:slider:4f15]
%   data: {"defaultValue":3,"label":"twist_root","max":8,"min":-8,"run":"Section","runOn":"ValueChanged","step":0.1}
%---
%[control:slider:32ef]
%   data: {"defaultValue":1,"label":"twist_break","max":8,"min":-8,"run":"Section","runOn":"ValueChanged","step":0.1}
%---
%[control:slider:bc80]
%   data: {"defaultValue":-1.5,"label":"twist_tip","max":8,"min":-8,"run":"Section","runOn":"ValueChanged","step":0.1}
%---
%[control:slider:e35f]
%   data: {"defaultValue":0.213,"label":"AuR(1)","max":0.5,"min":0.02,"run":"Section","runOn":"ValueChanged","step":0.001}
%---
%[control:slider:ac7f]
%   data: {"defaultValue":0.281,"label":"AuR(2)","max":0.6,"min":-0.2,"run":"Section","runOn":"ValueChanged","step":0.001}
%---
%[control:slider:43f6]
%   data: {"defaultValue":0.197,"label":"AuR(3)","max":0.6,"min":-0.2,"run":"Section","runOn":"ValueChanged","step":0.001}
%---
%[control:slider:5127]
%   data: {"defaultValue":0.268,"label":"AuR(4)","max":0.6,"min":-0.2,"run":"Section","runOn":"ValueChanged","step":0.001}
%---
%[control:slider:8529]
%   data: {"defaultValue":-0.101,"label":"AlR(2)","max":0.3,"min":-0.6,"run":"Section","runOn":"ValueChanged","step":0.001}
%---
%[control:slider:8ce2]
%   data: {"defaultValue":-0.156,"label":"AlR(3)","max":0.3,"min":-0.6,"run":"Section","runOn":"ValueChanged","step":0.001}
%---
%[control:slider:c1ea]
%   data: {"defaultValue":-0.09,"label":"AlR(4)","max":0.3,"min":-0.6,"run":"Section","runOn":"ValueChanged","step":0.001}
%---
%[control:slider:ecd7]
%   data: {"defaultValue":0.142,"label":"AuT(1)","max":0.5,"min":0.02,"run":"Section","runOn":"ValueChanged","step":0.001}
%---
%[control:slider:0ce0]
%   data: {"defaultValue":0.225,"label":"AuT(2)","max":0.6,"min":-0.2,"run":"Section","runOn":"ValueChanged","step":0.001}
%---
%[control:slider:f37a]
%   data: {"defaultValue":0.134,"label":"AuT(3)","max":0.6,"min":-0.2,"run":"Section","runOn":"ValueChanged","step":0.001}
%---
%[control:slider:0b54]
%   data: {"defaultValue":0.213,"label":"AuT(4)","max":0.6,"min":-0.2,"run":"Section","runOn":"ValueChanged","step":0.001}
%---
%[control:slider:83d6]
%   data: {"defaultValue":-0.029,"label":"AlT(2)","max":0.3,"min":-0.6,"run":"Section","runOn":"ValueChanged","step":0.001}
%---
%[control:slider:2f5c]
%   data: {"defaultValue":-0.101,"label":"AlT(3)","max":0.3,"min":-0.6,"run":"Section","runOn":"ValueChanged","step":0.001}
%---
%[control:slider:1c26]
%   data: {"defaultValue":-0.027,"label":"AlT(4)","max":0.3,"min":-0.6,"run":"Section","runOn":"ValueChanged","step":0.001}
%---
%[output:94a6dad1]
%   data: {"dataType":"text","outputData":{"text":"S = 71.91 m^2 (half wing) | AR = 9.01 | mac = 4.57 m | LE sweep 25.5 \/ 29.1 deg | dihedral 2.7 \/ 4.4 deg\n","truncated":false}}
%---
%[output:818cf68e]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAXoAAAEuCAYAAACAv9lxAAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQl4FMXWvRDIyp6ENQiIgCj4RET0PRUU1AeIPkXFXVweorg8kSjyUAEVN9wF9yeoiIrL764gAooiCggiixuyC4Q9JKBs\/3dqUkOnMz3Tk5me6e459f58\/Cbd1VWnqk\/fOnXr3ir79+\/fLyxEgAgQASLgWwSqkOh9O7bsGBEgAkRAIUCi50QgAkSACPgcARK9zweY3SMCRIAIkOg5B4gAESACPkeARO\/zAWb3iAARIAIkes4BIkAEiIDPESDR+3yA2T0iQASIAImec4AIEAEi4HMESPQ+H2B2jwgQASJAouccIAJEgAj4HAESvc8HmN0jAkSACJDoOQeIABEgAj5HgETv8wFm94gAESACJHrOASJABIiAzxEg0ft8gNk9IkAEiACJnnOACBABIuBzBEj0Ph9gdo8IEAEiQKLnHCACRIAI+BwBEr3PB5jdIwJEgAiQ6DkHiAARIAI+R4BE7\/MBZveIABEgAiR6zgEiQASIgM8RINH7fIDZvdRG4JffVkirls1SGwT2Xkj0nAREwKcILPt9tcz6doFs37ZVateuJW1aHyxVq1aVDn9r69Mes1tWCJDoOTeIgA8RWL9hk3w27RvZtm2rInf8pKWlSfXq6VK3bh1pfUgzadOKlr4Phz5kl0j0qTLS7GfKIFBSslPeeX+qFBcXSxXZHyD6tDSpllZNqqenS2ZmpmRl50jzgwrksDYFKYNLKneURJ\/Ko8+++w6BHSWlMuXzWbJp02bZ\/ddfUrVqmlRNC1jz1apVl4yMDMnMypbsrGxpe2hLKWhUx3cYsEMVESDRc1YQAZ8gAJKf\/d1CWbFitez6c1dQsgHZV6tWZs1nZEpWVrbk5edJi6Z5kp9HovfJ8IftBok+FUaZfUwJBBYu+kW+nbNAdu0ykjys+WpSvXp1SVeyTZaSbVod0kJaNK2XEriwk0KvG04CIuAHBEDy389fLNuLtxsseS3ZgOjTg7JNTk6OtDu0heTl1vRD19kHGwjQorcBEi8hAm5GAG6UX8+eL1u2bJa0qlWlSpmXTUCyCWjz6ekZahM2MytLGjdqJO3bchPWzWMa77aR6OONKOsjAglEAB42b783VYq3b5X9+wMeNpro9QYsZBtswmZkZkl2do5yreQmbAIHyQWPItG7YBDYBCJQGQSw+Tp56teyft0G2btvb0CyqVJFET1IPqDN4wcknyFZmdmSl5cvzZvWk9x6tSrzSN7jUQRI9B4dODY7tREAyU+dPluKijbKzp07y0geB6OqKJdKRfTwtFHaPH5gzWdLy5bNpVkTetqk2uwh0afaiLO\/nkcAJP\/Ntz\/IylVrpWRH8QG5pkqA5OE3j8NR1apXN2zCZknNGrWkTasCyatXw\/MYsAPRIUCijw4vXk0Eko7ADz\/+LPPgYVMW3iC4+VoFJ2APHI6CNq9OwmZgEzZbmhYUSJuW+UlvPxuQeARI9InHnE8kApVGAB42X3z1nWzbulWRetUqRi8bSDcHPG1A8hnp0OczJSc7Rw5uUcBN2Eoj7+0bSfTeHj+2PoUQQKAybL5u3bpZ9u+XgB5vInpswKpTsGXWfEY6TsJmSd16ueqAVL269J1PoSkT7CqJPhVHnX32HALQ5d9+9zMpLt4me\/fsLXOhrCJVINcYolMGiT49PeA7n5GhTsIe3PwgOahJbc\/1mw2ODwIk+vjgyFqIgGMIgOQRcnjDhiLZWVpabvPVeDgq4DdfZs3rk7CZWVK7Th1pcVC+5OfSpdKxQXJ5xSR6lw8Qm5faCGgPm2XLVkjpzpKgG2WVqmGs+eqw5tOVNo8AZk2bFsghzRjXJpVnEok+lUeffXc9AvCwmTV7npSWguTTlC6v5Jqyg1EqoQj85qvpA1IBl8p0hCPOyJScGjWkXp1sadv6INf3lQ10DgESvXPYsmYiEBMCIPm5835Um68geSXTaIIPEn3gcFQ1\/FSvXhbXRgcwy5Lc3HylzdetTd\/5mAbD4zeT6D0+gGy+PxH47fdV8tWsebKxqCi42RqUa4zWvAp1ENDmEbxMWfNlsg3i2jQ7qIkc1JibsP6cJfZ7RaK3jxWvJAIJQSDgYTNFNm\/aJPvLUgGWk2t0PBt1CvaA33wgSmWZt00WrPk8KWhYm+GIEzJq7n4Iid7d48PWpRgCIPlPpsyUtWv\/kL17dpfJNTgUdUCb1+6UIHntThnwtimz5stcKgsKmsjBTRnXJsWmUMjukug5C4iASxAAyc+avUB+W7ZcSksCHjYgeHUoyqDN62TfgU3YwAGpgDUfiDuPcMS1atWSxg1zeRLWJWOb7GaQ6JM9Anw+EShDAJuvX8+aIyXwsCk7CBUgeO1KGSD9Awm\/NckbfecDJ2Hz6zeQ\/LrpzAnL2aUQINFzIhABFyCwYOFPMmfej7Jly6ZyvvIVrflAdEodb\/7AJmz1ci6VTRo34CasC8bVLU0g0btlJNiOlEUAHjYzv5orRUXrDadetVxjOBiFMMRp8KU3etocsObhO5+VCWu+vmRV36OCmLEQAVr0nANEIMkIQJd\/853JsnnTRkMqQJC7SZsP5oG1sObhUpmRKdk5OdKkcSNpXkCXyiQPraseT4veVcPBxqQSAsrDZvKXsmbNWtm9Z3dZJEojwZcPcxA8BVuWPUrJNtWrS3qZt01mWVyb\/Nza0rQxvW1SaS5F6iuJPhJC\/DsRcAABkPyXX82RX39bLn\/u2mlyowxnzVczHZAKJBdB8m\/EtcnPry\/5datLbj1a9A4Mm2erJNF7dujYcC8jgM3Xr+BhUy4VYJkuX87TxhyG2Og7D5IPuFRmZmZKzZq1pEH9XG7CenliONR2Er1DwLJaImCFAEj+m2+\/D2SJUtp7QKIJuFFaWfNlCb+D4Q4CB6SQYERtwmZlS8OGDSRNdklLbsJy8pkQINFzShCBBCLw27JV8sVXc2TDhnWh3SjLrPmKmaM00R84IKVIviyuTY2cmtKwYX1p1oQx5xM4nJ55FIneM0PFhnodAejyk97+VDZtLJL9+\/dV0OWtrHmcgA2EO9DBy4ybsLDms6Rebp7UqpEhBzWp63WY2H4HECDROwAqqyQCZgRA8lOmzpIVK1fKHuVhUybXVKkiov5\/\/GOQbkK4U5aPUnlgEzYnp4bk5eVLnRoi+Xkkes6+igiQ6DkriIDDCOBA1NffzJfdu3fLtq1bAlo8\/qfYHf+nyT7wO2XZq9jzZRuxZmse4YhVlMp0gUtlndp1pF7dWswJ6\/A4erl6Er2XR49tdzUCsOKhyX+\/YIn8+ecule9Vb7hqQlebsCB39KSM5MuRvToFazwkhZjzgR8ckMrKxiZsI5E926XlwU1djQcblzwESPTJw55P9jECWqrZtr1Ytm3bIvv37S+LQHnAatdWvY5MGdToEbESkSu1RY90gWX6PKJVVleRKhGlMlNq1awluXn1pDk3YX08m2LvGok+dgxZAxEIImC04nft3BlI6K2tdpPFfsCqR4TKgGUfjDWvNHrEtdHWfFmSkbSyuPMZ6WUHpPIlM72qNCtg8m9OQ2sESPScHUQgTgiA5L+a9b38sa5Itm7ZLPv27wtJ8hUkmjIfeqNLZRqIvkybV1mkjElGINukZ0hOjZqSm5srtbL3cRM2TmPo12pI9H4dWfYroQj8umyl2nAtKdkhpaUl5Q8+BT1rDLKNwboPyDRGa768JV8uJHFZbBvINnXq1JUaORnSvICeNgkdbA8+jETvwUFjk92DgJJqflsp8xYskV27SqUUG65BEjeENDC4TwY3YNUp2IBsE\/SySTNLNsbYNlq2yZCcnBwV16bq\/hJp2qS+ewBhS1yJAInelcPCRnkBAZD8p1O+ku3FO2Trlk2yb3\/ZhqtZiy9zpbTS5I2brlqT1zHny23Clp2EzcrMlty8XKlWLU1aMByxF6ZK0ttIok\/6ELABXkRg\/g9LZd78xbJr104pLdkR3EwtR+YGgte6fOCglN58DcS5ObABe0CyCZ6ETaum8sJWL0v+DcmmZs2aUrduPfmztEhatWzmRfjY5gQjQKJPMOB8nLcR2LGjVKDHg+S3b98qe3bvLheQzOgyad50LX\/yNZQLZUC2CWjygSiVyqKH3zxcKjMyJFv5zTeQPX\/9KS0OoqeNt2dT4lpPok8c1nySxxFYt36jfDVrnqxbtz644WrU2zWRm0MaBH5\/QK8\/4COvrXmDJR90pzRo82UbsJlZWVKnTh3Zv2+f1M\/LkjzGnPf4jEpc80n0icOaT\/IoAges+EWyE77xkGpM5B0MZ6ClGbMurz1rtCtluRAHegPWaM3rAGZlp2AzM9UGbKOGDaS0ZLsc3Czfo2iy2clAgESfDNT5TM8gAJL\/ZMpM2V5cLJs3b1LWdDkd3uw6GWrjNajJV5RrtD4PicZ6AzYQoRKBy9LTq0nNzL1Su3ZNz2DIhiYfARJ98seALXApAvMXLJH5C3+SHcXbpSToG3\/AFz6owYfYdDV+DLQLZcC7xvr0q3kDFjlhkQ9WhTqoVVvy8vJk9cpf5OgObV2KGJvlVgRI9G4dGbYrKQhMnz5dPRd6PDxrEAIYh6AQ7x2EXW6DtRzBHwgzXEGTN8g1oV0pA2EOAidgtTYfSBOIE7DZOTnSsEF9KSktlVbNc5OCCx\/qbQRI9N4eP7Y+zgi0aNFCli9fHrLWuvUCJFuvXp7Uq5eryB+lVetDpVXrtsHQwwH3yQMbsBU3Xw8ciirvZWNILqJ95rOy1QlY6PM7i9fLIS0PinOPWV0qIECiT4VRZh8tEQCpjxgxQv09v34DefCB+yuNFog\/NzcvQPpVqkjrNofJoYe2U6SvT77qU7ABTd7sTplW5jMfsOYzMzKlRs2a0qB+vmzYUCTtDm1c6bbxxtRGgESf2uOfkr0HuY8fP171vWjjFhnz5GMx4wCS37xpY4V6cvPypc2hh0teXn05q88F5aJTmjdgcTAKuryKNZ+eIVnZCHOQJ\/v3i9TM2if16nIDNuaBStEKSPQpOvCp1m1o7zNmzJCS0p2yePFS+fCDd8NCAIJGbleQ9KaNG9S1uXkHYsr8tHRRpSDMy68vJ3bpLuf2vSQQnVIfkALJpyEXbDUVaz4zM1Nq164jdevWlR3bN0rL5g0q9TzeRASAAIme88C3CECSadasmZSW7pJx48fJd9\/Ojkju6qUoI3WQPAjdTOqwzv9xwknqQ3Bo28Nl48YiWbrkR\/WjCwgd10HCWbJ4YYXnQibqevJpcsGFlyu5BmGIq6kUgQHJJqdGDalfP19lpcqtU01q1azh23Fix5xHgETvPMZ8QoIQgCSjvWZyatSR8849KyKxg4xhscNah586SNtI7Ph7Xn6+5OU1kDZttXVfpRyxg9RPOLGb+gGxfzHjM\/nyi6mysWiD4G+HHXaE+heEbyb9+g0aSvdTesql\/a5WoQ4Ckk22imUD3\/mtm\/+Qtq25AZugKeTbx5DofTu0qdExTezo7YpVf0i\/Sy+07DjkGBA3iDygm4cmdljpsOs1sW8sKpKZX34erPfQtu2k7WHtFbFv3LhBkfc7b01Uf9fSTJ9zL5KiovXy5Yyp8takCepv+F2Xk06RGdM+k7cmvRKsr0GDRspt87Zhd8uxx50gtWrVUj7z+FC0bdUwNQaSvXQUARK9o\/CycicQALk3b95ccmrUluHDR8jYMdabqZrctQwDuQXkHJRYDDKMJvalSxYFZZiAtX6yHNq2veTnN1CWOiQabZmffc6FcmLX7upvRmLHf+P355x3sfpofDF9irz5xitKF+rS9VQ57\/xLZNLrL6sfYxl5z8Ny5b+vlb1798mfJUV0p3RiAqVgnST6FBx0L3UZcoz2a+90TGelt\/fq1cNSbzcTu1mGwd+RsKNN23bKYsYZqFDEfsKJ3RVMIHajDIONVG2tfzF9qnw54zNF8G0PO0KuuW6QInzUO6OM2KHFX3vdYDm83d+UrDN92mQZ8\/iD6r\/vvvdRKSraINOmfiITJ7yonocPx8ndTlEbtu3bFnhpqNhWFyNAonfx4KRy04ySzEknnRQRin8c3zWkvq6Iv2xjFJVgo\/X\/3nn9gEVfZrGD2EHYsNaNMgyIt0vXU9TfQOxvvzlBkbnRWt9YtF5Z6yD3c867RFnreBJ+D2J\/47WX5Pobb5Vup\/RQces3btggjz0ySrXh0SeeVy6VPy6cL3ePGCpr164WSENfzpzJ6JQRR50X2EWARG8XKV7nKALjxo2TFStWKC+Z+vUbyvARwyN6yZgbpL1hIM0cf8JJylKHRa+9YUD4x59wstLXtQyDa6Gja2v6hC7d1PMVsc+YKm9PChD7gIGD5LDDj1AEDn19xrQp6qthtNZxz\/Bhg1Vdz\/7vtUCKwCpVZfGiBXLr4IHy6hsfStOmzSQjI13504+481aZ+91sef2tj6W4uFi2bt0iF57XSxo1LpBXJ7wsXbt2dRRzVp46CJDoU2esXdNT+LKf3qunak+nTp2lpHSXTJr0mhRtWB+2jbB0jS6MkTqkiV17wxj1db1pevY5F8nGjbDWP1OEH06GAdEPv2u0kmCg+Ux67aWAtf6fIcpzJg35XtPS5LUJL8r69X\/IsDvvlYz0dGWxr1mzSnr36CKjH3lamhQcJJs3bZJVq1bIffcMkx69\/iXYkF254neZ\/c1MOffcc6RLly7Sr1+\/SF3k34mALQRI9LZg4kXRIKB1dWyYQqOGrl4\/v57aQEVpfnBrmf755LBVGuUWs3+6tshRQSiXxVAVoz67MsxTYx6SxYt+UIeazu17qSJ2yDBPPv6A+hjdc99j0qhxExWIDJu7ox8YKb3P6CNnnnWuuvavv\/6SO28vlDWrV8m11w+WLZs2S9HG9fL2m6\/Kvr375MijOilSX\/77r6qpNWrWkjatD1F7ESB3rGpI8tHMOF4bCQESfSSE+PeICOhTp9pCl6rp0qtHt+B9IPtI1jou1m6LRq8W\/B7Eri1wELu2vLVWDr0cBVY5fmCVm4v2ijHKMCBlaOv4QRuHjxwt8GvH70H0d\/x3kFxw0eVy0SVXBkMTYNN07pzZctc9Dwn0f0gu\/\/f26zL2yYdkUOEdsmfPbhXWeNbXX0idOrWluHiHIv6ff1qspJgXX3xREbreYNaEDgzxO1yjP4gRgecFRMAmAiR6m0DxsgMIgJRw6lRZ7PkNAm6Odw6NGiLo4dDVUYzkrkkdBA7ihhcLDiGhwOPlxJNOUZuq0Mk1sUNWuWbgzYqwizZskDcnvRzQ0S0KrPXzzg9Y6\/iBBDPt80\/lxZfeUrFmEIZg8+aN8vyzT0iHo46RriedouSWLVs2y223XC+NmzSVTsccp+SXRQvny9GdjpauXbrIjz8uVIQNbEDqLETADQiQ6N0wCi5rg7Y4YV1CTx879mm5pXCQIi8QKqQYeJNEKmZNHQQO6xw\/mrzxLz4WIP0+512sJBKQt\/mQEXY+4dXyxbQpymKHfzo8XCCXY7UwCV4v0yaXyS0Br5eionUy9onRsujHBWUfgIClj83WDRvWKUv9iqsGSmZmhpJhPvzgHXn+2Sdl2B33yu7de2TVquXy6cfvSxXZK+ece4GMH\/c\/2blrp5JZnnnmOZXtSXsH0RKPNBv492QiQKJPJvouebaWXrAB2LxFaxUXZsSd\/1Wtg4UMIoaUEamAyHGt8Zg\/dPETunRX\/uZaT9eSS8AffYM6TKQt8z7nXizn9L04aLFDVtHEDis8QODrZewTD6p\/B15fKIe3OzJA+EUb5M7\/DpJ27Y+UG28eqjxeqlVLU\/sBH3\/0rvQ8\/Sz59OP35Pt530n7IzrIuj\/WSIeOnVUo4TVrVsry33+TmjUypbDwVlm9erVMnBg40Tpt2rRg17W1HgkL\/p0IuAkBEr2bRiNBbYHsogmrWfNWMn3GdBn\/4nNRPV2TOqxzbDKi5MG\/vEs35WMOSx3EDsvc6HcOcjZKLgHLPHB6dPHiH+SpJ0erurTbovr9ogUy5gkcMjpSbvjPrUG3RcgtPyyYK0OGjlSbo\/v27ZM9e\/bIA\/cNl9KSEul7wWXKLx06+vvvTpL9+\/dLTnaO\/PLLUpW16ZhOnWiZRzXqvNirCJDovTpyYdoNC137YJeU7JTLL79cvvtutnJlzK5RS37\/\/XdlQUdbYJ2DqDWJByQXHPMPxHXRkov5QFHQ71xErlGnROGPLmqVMPbJ0XL44X+TgTfcorRy5FSFDDN1ykfKFbFx46Yisl+R+KAb+6tEHJf1GyCrVi5XLoufTflIOnQ4Uv7Z40x55+035IcF8yUnJ0tGj35I2rc\/ImiVGzdAgQ09W6IdfV7vZQRI9F4evbK247AR9HMQ+YaNm6XvuX0UsUN2walOaNtWXi+wzFFA3sYCzRyErT1mYLUfkFwuVr7nsMy1xT7gupvlcHWgSIKeLNDz74TfOez1Mg+X1yeOlxfGT1IEjkND0MZvuuEqaXlIa+nZ6yyli69Zs1pmz\/pSWrRoJmecebYs+nGhvP7aRFm3bq3MmTNXatWuI9WrVVXJQ0Dgxk1P7c1CzxUfTGx2IW4IkOjjBqVzFWki11a6PkXao2dv2V8lXXr3DBzfRwExg2DXr\/vD8nARSByyiiZ3vRl6YtdTZOnihUpCgXVudEVUkgs2Sd94RdV\/x4gHlS6uiR3Bud79cIbyWAF5r1+\/Th56YKS6tvupPWXVyhXKr\/zrr6ZLp6M7yHF\/7yIZGdXl448+lNmzZykLHP0rLt6uNjhxSvbOO+90DlTWTARSCAESvQcGW53ELCuwVI3Jq0HGhx7WTlnsoU6Nmr1c1Mcgv4FcPfCmsoiL2AyF\/\/kU9d+Pjx0f9GTRgblG3P2QHHV0Z0mrWlU959UJL6p46edfeKls37ZVpkz+RObO+UaO6XSU9Ox1puTm1pWlS5bIgw\/erwi8T5+zVet1vHh6qHhg0rGJvkKARO\/y4YR1ayeoF7oBGQaeLAEvGRwsChwewn\/fPvz+oAwDAofkguP87Y84SrIy02VDUZG8PnGc7Nq1S66\/YbAUb98uU6Z8Igvmz5HcerXl9tuHy959e+WXn3+WSZNeV9a2XmEYT8JSMnH5hGLzUhIBEr3Lh10nsh4+fHjIlsJiB7lDisHmJv4dMPBmJZdkpFeTdevWK7klLa2KnH5GHylo0kTWrl0jo+4ZLu3bt5OB1w6U3Lx6yhXxgw\/eVzFWjMG0qHm7fIKweUTABgIkehsgueUSuEVCnzdKN7Cg7x71kLRocbCUlBRL\/6sul65dT5AzzvyXHNziYKldu6batDR7mZDA3TKqbAcRcB4BEr3zGDvyBLhMgvThccIAWI5AzEqJgG8QINH7ZijZESJABIhAaARI9JwZRIAIEAGfI0Ci9\/kAs3tEgAgQARI95wARIAJEwOcIkOh9PsDsHhEgAkSARM85QASIABHwOQIkep8PMLtHBIgAESDRcw4QASJABHyOAIne5wPM7hEBIkAESPScA0SACBABnyNAovf5ALN7RIAIEAESPecAESACRMDnCJDofT7A7B4RIAJEgETPOUAEiAAR8DkCJHqfDzC7RwSIABEg0XMOEAEiQAR8jgCJ3ucDzO4RASJABEj0nANEgAgQAZ8jQKL3+QCze0SACBABEj3nABEgAkTA5wiQ6H0+wOweESACRIBEzzlABIgAEfA5AiR6nw8wu0cEiAARINFzDhABIkAEfI4Aid7nA8zuEQEiQARI9JwDRIAIEAGfI0Ci9\/kAs3tEgAgQARI95wARIAJEwOcIkOh9PsDsHhEgAkSARM85QASIABHwOQIkep8PMLtHBIgAESDRcw4QASJABHyOAIne5wPM7hEBIkAESPScA0SACBABnyNAovf5ALN7RIAIEAESPecAESACRMDnCJDofT7A7B4RIAJEgETPOUAEiAAR8DkCJHqfDzC7RwSIABEg0XMOEAEiQAR8jgCJ3ucDzO4RASJABEj0nANEgAgQAZ8jQKL3+QCze0SACBABEj3nABEgAkTA5wiQ6H0+wOweESACRIBEzzlABIgAEfA5AiR6nw8wu0cEiAARINFzDhABIkAEfI4Aid7nA8zuOY\/A\/v375YsvvpDJkyfLEUccIV9\/\/bVcddVV0r59e1sPLy4ulptvvlm6desmffv2tXUPLyIC0SBAoo8GLV5LBEIg8NFHH8nDDz8szz\/\/vDRv3lzmzp0rt9xyi4wdO1batGkTFjN8JF566SUZMWKE3HvvvSR6zjBHEEhpov\/mm2\/UC8nibgRuvPFG6dOnjysbuWXLFhk4cKC0atVKbr\/9dqlWrZps375drrvuOmnYsKHcddddkpGRYdn2H3\/8UV2zbNkyGTx4cESi55x15TSo0Ci3zdmUJvq33npLCgsLpaCgwNbsWb16dVTXotJo6s7du1e6ichPBQWyKS3Nsk1oR7R1R3t9NO12uu4HH3zQtUT\/3XffyUUXXSQjR46U888\/X40LrHRY559\/\/rm8+OKL0rRp05BjCclm9OjRcvLJJ8t\/\/\/tfuf766yMSfTRz1u58tTuf7F6Hzsb72U7UGW1\/7MxzXafb5iyJvrBQ7A4KXjIUu9Ylru\/cubMtsoeltuaTT6TP3Lkio0aJRNB3UXeTJk3k2GOPjfiReuONN+Sxxx6T119\/3VZboqkbE3v27Nm2+xkNhprU7I5PRCAcuACY3nbbbfLEE09Ir169gk946KGH5JlnnpEJEyZIp06dKjwZH4OXX35Z8vLy1Grg8ssvj4roX3311YhjbxdrNffWrLE1r+3ODbt1Yj8DxtagQYMiPt9uf\/SctPOexrs\/6PeFF15om1McmJIhqyTRFxaKnZcmIQOycKHI0KGhiX7jRpFvvxXp2bN8U\/D7G28Ueewxkby8kM3csWOHlJSUKFJJC7NSSEgfo3iIW18aYxdee+01GTp0aIU5pD8A0N+PP\/74Cr2GZPPuu+8quWblypWOEH0UUCft0p07dyqpq1atWpKVlZW0dsTrwW6dsyR6LxH96aeL9O8vcsUVIqUlvAKVAAAgAElEQVSlIv\/+t8irr4rceqvIffdZzlUSfbxe44r1VIboIdncd9998u9\/\/1tt3v7yyy9REz1WildffbUcddRRznUuATWD6IFHzZo1PU\/08+bNkw8\/\/FCwSnDbKpRE7xWix0tnJHf8dwSC1+8pid45xrKSbp588knliWNeLRolm55lq7PKED16NGzYsHJykXO9dK7mXbt2BYk+MzPTuQcloGaQ\/N13362eRKJ3APA9e\/bIe++9pyyj\/Px8Wbp0qWDXGxp2uKI1YGh5GJikl3DSjW6cJnts8IWx4o198SrRQ7t1o3VkxPb777+XSy+9VG2m2tmMhQU7ZMgQef\/990NOt8aNG6sNXOj2oYqes66RG2N8aSjdxAigzds9b9HDQnr66aflyy+\/lDFjxkjdunXVSwRtFBtkcHGzKvqlwTJ44sSJNiFz8DIz0UN\/h1wze3b4h3buLPLBB77T6C+44AK10es268g4GKHcK\/FhxeZinTp1KrhXYr5CqoBxosvvv\/8uN9xwg\/Tv31\/OOOMMJWPATZNE7+C75lDV1OgdAnb58uVy5ZVXKhe3K6Bdi8iqVauU5nnmmWcqf+YqVap446WxY9FXAkevWvRufWnMQ4ADU5BqcEBKH5iCb\/2jjz6qPGPg0QJf+RNPPFHw8TLPx8pIN7ToK\/EiJOAWt85Zz1v077zzjjo+bvRugIcJlsfbtm1TLyB29D1hHZHoyw2TW18a81zat2+ffPrppzJz5kzp2LGjzJgxQy6++GI5+uijFanDYr\/22msFmvw111xTwVqvDNFDmsSP14vfpBu4MePHbatQzxP9\/fffL6+88ooi+g4dOqh5\/+effyqXtzlz5si4ceOkRYsW3iD6SG\/tjBkiXbtWvMqn0o1XiD7SsMXz766TG2PsnN+I3q1yo6eJXp9AxNLZvIGFDwA2aPH71q1be5\/otV5\/xx0BPf7660WaNQu4WHbrFnC5tCiUbmJkIxfdHuoQGcgF8pAXy969e1WzvXC+A84dkfby3HrIj0RfWKg2xrCshl+y1XH1ZL9EVTdvlhoDB0rJiBGS8dprsu+QQ2TXhRdK2s8\/S86dd8qOMWNkX716IZsJKau0tFTq1avniRcKnYBPMlzVcMrRbcvgZM6FUF430P5ZnEVAh3RAlNJwxa2rUE8TPQAPJd389ddfcscddyjN1I6rGnx5oefjX+j9ONLuulJaKjn\/+Y\/sa9JEdnfvLjV795Ydr78u6ZMmSdXff5cdr70m+3NzLYkeZJ+bm+spokdERxJ9+SENR\/SRSMh1c9pDDdIf00gYk+gdGlQcUkAwqFg2Y\/\/44w9F8sbStWtX9ZGAF4WrCuLgnHWWyIYNAb3+wgtFnntOJDub0o2rBsqZxoSSBuySkDMtSo1a7WJM6cah+RDKvRLEDZfLHj162HKvRMxwuMFh4xb1GQuIHoQP4ne80OumHMRutY4cnwdhHhBqM9YuCSWz3V5\/tl2MuRnr0EjrA1OzZs1Sbk36wBR8mJ999llp2bKl5ZP1SwOi\/\/jjj9V106dPVz74oQj\/zjvvlH79+jnUExEh0ZPoI8wuPWeN7pV2Sci5iev\/mu1iTPdKB+cC3CknTZqkPA+aNWum3CrhrxyO5NGcUESvmwmiB+GD+M0WPsgepB\/3QqIn0dskeuOBKbskFPf5mkIV2sXYratQz2\/GxjLXQi2D4XNvJHJN9qEs\/LgTfiSiNwY1gzZ\/000i+OCMH28Z\/gD4eNW90q3L4FjmXKz3cjM2egR\/\/vlnmTp1qjL+KltI9JVFzgX36Zfmb3\/7m+CELYqZ6I0W\/vjx42X48OHOWfjhiF6TPHzmIUc9+2xgE\/a110SmTg27IetVonfrMtg8dSuTHByrT8Ri0pLh6aefrpwKwsVmMq5CadHbI5ANGzaoOEI4THkror1WspDoKwmcG24LJ91YtQ8btrDyzQWbtjFb+JESj1x2mcjo0QGPG030K1aIDB4c1qr3KtG7dRlsHvtok4MXFRWp\/LKIxdSuXTt588035amnnlIb\/shMhaBmVsWvm7E4IYvAhC+88IKcdtppgjj\/NWrUUAELEUYCp9\/r168vv\/76qwohgQBwiFSLAjxD\/R1B5ZBd7fHHH1eZ1ZAABpgffPDBUdOPXaJ36yqU0k1hoRg3Y80zQEs3iFdiLJBy4m7hR5JuhgxBxDaRiy4SmTAhIN2cc44IcpUy8UjUL288bqhMcnC4BOPj27dvX9UEHNhDzlmsKo2hPEK1z6+bsSB65NjFquZf\/\/qXImSQerdu3VQ4ExBo7969VXgTfCRxvgIrPhRsTFv9HSdu4YEHr7pEWPRuXYWS6CtJ9EZJ56STTqrgpYO\/Q+a57LLL7PviRyJ6VGqOd\/PCC2HDH+AWWvTxoPTQdVQmOTiCnrVv316dVNYFCUxAYFY5ZvV1djX6I++e5VynY6z5Hy3ryJgL2laoRa\/gjGdi8FEcNWpUuVAm+CBcddVVSvpCCff34447LqFE79ZVKIm+sFDsaPR6VsILB8QOC994mMrKwo+K8O0QfSVeMq8SvVutI+MQVDY5uHkYIVWA4CBdNGrUyHKU7RJ9vUHTKjFTEnPLBZ0ahiV64\/4DpCysdIwn3LG5ipX22WefrRoc7u+4LpEWPYk+MXMoqqdUJhKgFdEbH4wPgdktE3+PeNrWrkbftqI1FK7jXiV6t+qdRqwrkzPWPFZajoAjwIABAyzzJ+A+u0TvZYveSPQIMw6JFESP\/QwUHdZZZ\/QK9\/cLL7yQRC8itOgLC+Xwww+3TO1mpdEbX1YrTx3EaonraVto9IjEGSZSZSjC9yrRu\/U4ebyJ\/scff1R68z333KM2HMMVjQk2F2Gp9urVS+nZ0KIjxWGJygpKwsXaGjYSPWQuJFF\/5JFHgvlxgRf6\/sADD6hWhvs7JLJ4WfSItPnuu+9aIuPmQHwkepNGb+VHb96MtUP0+ppwnjq2T9uGSyvIePRJoKXAI6NNDm5uKD7CiM6JDGlW4bSN9xhXodCpjzrqKIFrpl+JHmkXEWgQRXskQeJasGCBwJCCa2u4v+\/evVt56OC9Rl5fbJ5XJn6V9rqJRPRuDcRHoo9A9OYX0450oy35UJ46CT9ty81YRz8C0SYHNzYG3jZwHzzyyCOV65+dYle6sVOXm67ZvHmzygYHowhuyjjcpN0nQc7YeIVFjdPucKdE2lCEO0EJ93d8CP73v\/+ptI5\/\/\/vfVVRbxJWPtth1r6RGHy2yCbi+Mn70sRC97hI2bkN56sTFFz8Ebl6Vbtz60hghjjY5uL4XJA9iA+mcc845QV0eboNIfRlN+ku7JJSAV8q3j7CLsVvnLC16k9eNeabGotHruqL11LEkfGMIBF25j8MUe8HrBsMQbXJwWJmIzYSPRJ8+fYI5ArAChDvhoEGD1GGhUMWvFr3bvxAkerePUJj22fG6cZLojRa+laQzbdq0gKaoSb5p0wOHo\/TvUFGYmPRetei94HUD6KNNDo4Ul0heD6veXPBxw8Egq+LXA1NupxG7RO9W44QWvcnrxsqDJpJ1bpyoVhq98ZpQz1n+\/vvyW9++ctPOnbLQcDHcMl96+GFpOmxYIASC0b1yyRLfhkDwgtdNognKryEQEo1jtM+zS\/RuNU5I9BFOxponRDw0etQZ8oNS5ke\/vH9\/ufzhhyv44v\/fGWfImevXB5KD5+UdsPL79xfp0sVy7nrVoner3hktScTzeko38UTTfl12id6tc5ZEH2MIBEvSNswhux+Hj+6\/X2YjOiXSBbZvL6u+\/152dO0qbbdvDz8j6V5p\/431+JUk+uQMIIk+ObhHfOr27dsFPrjhXKlCed044UcfDdGPEpH5ffvKzS++GOxjqKxX0XjopLJFj0TxSIxeu3ZtqVq1asR54\/YLSPTJGSESfXJwr\/DU999\/X0Wx0wU+uPDJbRsmXIB+aTp27Ki8IOJpnePZUR+yWrhQpnburIi+42WXVYipg49QZXLaepXoCwsL1ZF\/HCiCd0q0BV4tmBMzZ86U559\/Xk4++eRoq3Dd9dyMTc6Q2CV6bsY6OD4gMoR5xaETHREQsagPO+ywsFacHa+bZGj0kG6mb9oUVfC0cHF0vEr0sW5sIbIkjssjDPUhhxwSDC+AkLhjx46V7t27q4B2Xip+3ozFPMUBNIR3gFHjZMGKH6drv\/3223KHr6yeaZfoY52zTvXZFxo94mHgB3GrkZTAbrFzYCoR7pVoLzx1zBp9uH5Y5bRFaGRzPluvEn2sG1sg+j\/++EPOOOOMClDC2h89erT85z\/\/CZ7AtDtvknldKE8kuySUzHZHejbGAwlCEKDMGOsm0n2V+TvOMWzcuFEZhzighoihxpDRoeq0i7FbPcU8T\/Q68h+y9ECPhw\/yxRdfLI0bN444B9xI9KE0et2RkC6ZFknMjYSfqkS\/bt06dXR+8ODBwePyxkmByJOZmZkq0YVXip81+lBBzZwcl\/vvv1\/wzHgSfazGiVP99TzR6+Xe119\/LZMnT1a6OAgRqdkiBYmysxmbaOlGa\/TGzdhwRK\/\/Fm7DFkGfsCGZl5cXPIXp1ISKZ72xvjSw3BBLZvHixSqglTkfK4geS3gEvfJKIdHHb6RI9PHDskJNWDKhgHTiXWDdY1P1rrvuks6dO6vQrzrwUahn2bHoE030MnSocq+00uiN7Ql1MAu\/MycwP+iggxQWWO0gyqFXSqxEj37CEICU9cEHH6hkFQjp26BBA1m4cKFauuvUdYnGBKdicUIWsdXhOLB06VK1cRwp4JZfpRvgr8d74MCB6gwJEoz07NlTSbIICWGVU7Z69eoqcBnkH1yDrG6XXnppMKInDEA4XLz99tsq5SDqRB5aI9HDOwupBrE\/gHtbtWpVbn+P0o3pDUF2eyQ3DhWUCZMbS2mk90IORyeKtuIw8MZkBeGIHn8bNmxYMN618dqbbrpJZs2apSahVUFfzj333GC4VPN1uB+Bq1BHU4QwCFEQgvXLsWNl+t\/\/LjuGDpU9NpKL4B582MxtW7VqlUqK\/PDDD6snYVLDtRAfP\/zOqg1OjEdl60TMl7vvvlvdXlmvG\/1sGABIHg0NGC63uuDlR1TDSORa2T5Y3afnKEhpzJgxyhiB1xjC7+L9MK88jPUYN2PvH\/OiLPtpqfT7d385qF5muXj0xdPHxbvZcauvWn5zyTq8a4X6NNEjeiX2Tn766SfB+4f367bbblPz3JxTFiGaEau+R48e0qVLF4UjVrDPPfecSgKOjzvCOSNJiTmLlyZ6XIsPy65du9R1ofb57BJ9ymzGIj43ynnnnVcuUw5Ci+KlmjhxouObLfC0wADjpQn3QdEvTW5uriITxPZGooJLLrlEbrnlFtWPa6+9VhE9vAGsivke83V4obEh+MMPP1iS7L333qvIGR\/DE044QXmFmIvd54AkEZIVhA\/rBR9fED3CvKIgYQOSNbi5IInDM888I7Nnz7ZF9CBPyFPZ2dmWnlaw2pYsWaLIHsYIXG\/T09MTDgM20pEMAzHoEW4XBWOlVxzXXXedZZYpPWcHXnimbN1bS65cN0ouWthU0rILyhH9snOrJLxfdh9Ys2s\/yR944JyIvs+s0WNMEeETnAKjbdOmTYKMUcacsngvkTMW7wtWRvrdh2GFlQGseLzXyE4Fokd9Oi0hiB6S72mnnaYMUyuSR\/vsEn3KuFdiUO677z71VYWVi2UVPGKwTAbh5OTkKP3cKYseg4LlMDbgQN7hdHrXafRloRHgKonlpzk3rd04PAiEhjp0gSwAYoc\/uS7RHLiy+wLH+7popBu4TCJQGF50SFSQaJo1a+bKQ1LIcYr3w0hY+Eih\/du2bVNkFClM8UOt\/5B\/9Kkj1RpkSu+RVXxJ9JhPOhE4PHFQQPRGrxzwDTZT4SoLrtEFq1YYV\/gdJCB488DgwsfDSPQwJmDBg7N0DtpQ8xhEv3vDcrn8qrPlrJN7S92mLQUrE3OJZs7G+30JV1\/cN2Phm4pTiCAVpDbDhIUWieUxvqCwniElYFnlVMGHBYPvOY3eAIid07R28dNeNyBBLGONh67cTPjRvjTQaGGhQaaZO3eusvBgbGBZb9Zc7WLnxHV4D9BGEH2HDh3UIyAvQYueM2eOOuhn5UeujZNnrtsvHQ8JtC4U0XtZujESOeYAjBRgBYzMRA+VQL\/vZhdJkDq88fBhxSr2q6++Uqt8I9HjPcPKDrwVzoEDRL+3dLW8f8d+hfm+HXtk9a48abhlj\/x56q2yv+6x0rzgkOA+Q6xyY7znnSNEjxcM+hqs+BUrVqiDTJFIt7Id+\/jjj1UeR2SkOeKII9SyHBuSZ511VsSsPXY2YxPpRx\/pNK0ZI7tRMrGMxQpHe92E2rAF4eMFMK4EKjsm8bovWqI3PhfzAIaGJn0YHyB9WPsg\/WjOW8SrP6gH5AOZDjHsNeHo+vEBgFGE31utRM1En9H0HOlx42xVhR9zxkKGfPbZZ5WMBwPFTPQgcUg3OPmsD78BY2CBvQ9o\/Jjv4CBY\/2aixxzD\/bgGvIWxCZULwEz0VnPi+7UF8u8H1tiSG+M5ryLVFXeiB5g4iAI9rVOnTnL77bfL\/PnzlaUFFzedHixe2igOxeCDguXZP\/\/5T8GJWLzQSPiNlUO4ol8ap5KD62fbsc71galo\/egjDTAsQzPR4x68NFjOmj10wp2wjfSseP89XodPQPqw5rBxjeV7Mkk\/XkQPrIcOvlJ697lKSQ5+yhmrJS1Y8Jif0Ngxh7FfYyZ6zGO4x2JM8aFs1KiRCptRWlqqVkvYB4EXHs5KwPqHhY+PAk7NwzFB+9HDBReaPpK+wOvGzB1KutmxWsbdFLDoQ5Up80UefifAOb636EH0yMsIDRI6c0ZGhuo4iBhfS1hSOit7vIkh2vpCHSe3q4PDX98qybBdSxveBeVOsRpi3UTrR4++6w+KWaM3HpiCJ4tx5WCVuBwvjd6QjhbXeF3vhAcDsMCx9wkTJijSx1yEW10iSyjpBhvFeG8geZotfWPb9Jw9\/RiR045rKMd0v1r+1f9ZXxA9pDdsqoJ8YYHjwBscE6C147yDVU5ZuMqCWLVXGwhb3wNMsXrCniD2bWCEIhzGgAEDVH04OX3PPfcowxDSGWQgOGBgJW\/cJ9EW\/auF1kT\/w+8iL31eReb+mgJEjy8mADrllFMqfBUxWPhCwyPHyc1Yuy+tfmmw5MMS0E6xa51HkmHCxaOPxY8+VB\/CEb2+PlRIhWTr97F6MOBjjCU8rEE4AZjLypUrlScUNr1D\/d3OfKjMNZAj4CYYy2ZsJI2+Mu3iPdYIgOhXbt4lb15ULFX2NZS11RpI86Y1pGp2gdqUhXvrH9Uayup5c6XfbWP8b9GH86MHjE4emIp2orpNo48m1g36anflEEq6scJKf8iMfwfhY+WBFUgiSywaPdoJ7ytYZkceeaSy1sxhMWDpQfvFS4zDSonS7UO5V8KyxOoCG8d23CvHnLtTOv8jUw1HqM3YRI5TKjzLrntlrHPWKSzjrtE71VAn6rXjXpnIzVgQfSwafSjpxkqj13hayT1u2LCN9aXRRA\/iXLZsmdqgN4etxioTB2auvvpqy4Tc8Z57+sAUzmdoJwUc9IGGjI3Hli1bWj5Sz9l6Z9whj66\/TQ7K3SoXfNOogntlvNuc6vWR6D08A+wQvbl7Tks3keLRV9D1DQ20o9GbQyBY3YNqk71hGy+ih04LrRcut9iEw4YzDpBBE8aeBT6G0HUjbd7Hc6rrcB1YAcPfH26V8BwLR\/J4vp6ziM9\/5tX\/lXtemSZL3ry3wsnYeLaVddk\/MBVrDgWnsKZFX1ioXnh460CigIum9sfFC6g3XPEvfpwm+mhi3didFHY0emNd5v0DEL6WiYzXhQqJbLdNdq6LVaPXFj2I\/tRTT5V9+\/apk5I46o7TkHAYANHCC8Pp+Od2+mvnGj\/Ho7fT\/2RdY9eid8KBIB59JtEXFqoXHxqp3WImf2y8wkrETywfBzdo9MDAyvMomhj4drEMd12sL83atWvVBwpBzLD5CddeyCZw00PcG1jRN9xwQ8Kt+Viw8XP0ylhwcfpeu0Qf6yrUqX6Q6AsLlWsVTo2a0\/TFG3S9KsC\/CMCEw2Qo+DjArWvuuHHySFZWhZyxuh2RXD9DtbeyGr2xLvOmbyiXTCcOXMXLjx6ui\/gxHoRBbBn4TOPUJax9rxQSfXJGikSfHNzj8tRQLw0q1oSPf\/GjCdn8+3h\/GNqLqM3YZ5s3l479+qmPgfGkaiSit6PRm\/3ow2n0GmQr7x6nN2ydto6wEbt+\/Xp1UtYrhUSfnJEi0ScH97g81Yroo63c+AHQHwrjx8H8gbCqXxP9UBFZaLhIy0LmFIF23CtRTbQafTT9d3LD1mmij6afbrlWz1m4g+IHxS4JuaUPXmyHXYxj3VdyChtKN4WFjodNDjV4VqsG6Mg45We1WtCHmHDqGCEM7BzMisaPHm2NtHJA23CN8QSuExu2JPqKMyeVNmMht+E0PfZSkCks0XkDjOjbJfpY95VI9A4gEMo6cuAxlapSW8qQVvATqoD04XuNOB+w+uPpRx+u0aGIXl8fzxO2brWOKjWgcbrJz9KNjjaJjGhIkoP\/xlkD7KcgpEFWVlacUIy+GrtE71bjhBZ9YaGaVEiI4tZil\/RB9rDww8W6icaPXuNhRyIyrwKs9PtoTti61TpK5jzxM9Ejzg\/OM+C0shtCpFTGoifRJ\/PtsHh2vDT6RHZNbxDDS8ccedLYDqOuH61GH0m6CdVfq3ti2bB160uTyPE2P8uvRI8MdAghjLDSiMaJg22IiZWZmalCWUC2gZHy6aefqnDDOACHiLiYI\/gowMnAyRwXtOiTOetjfLYXid7cZW3tw5K30vWxFO7bt6865h\/p5GVlNfpIUk+oA1eRQiJ7heh1\/HMkocaBOxAQErwg9aNVwUEtEBvyKaAg9yl8\/cPli8V1djdjt352Qoxvh3O3V889VnI6PFjhAeZUgohKCY0ewed0RE9gjMiTwBeb0fgbQlzAcwqRKUPFko9HT+wSvVvlRko3SdqMDTn5Fi4UfTJWwpCE1cS1knh0cnDkjIWuj8mICKNWur65frN0E06jDyf3RLth6xWiR3A0xDbHCVvgi+xWCPGMkLtt2rSpMFywXpGnAWFzkcsUEV+R3QjjgYTvyGdrVexuxm5+r0U8+M2ROpAsxQ7R4+HmZCFWeWWRDzZcwpZYO2KX6N0qN5LofUT0xslsJH3IPMbk4Po6kJKVrm+sy45GH+njYG6bDhZnbIs5jo8XiB7hM6Arw6IEeSMCJvzzYWXCOkfSC52TQfcVYYohp2GVhYKE8Eieg1DZxvSCocgp1CGyUCTkB4veDtHjGmj7SBZiTEEYK7Gb77dL9PE65Bfv9nuG6JE8GckiEIUQG3rGHXi8KEjBBi0P7olIho1lXSR3LGOAKLg0Jr3EaNFbtR944GfevHmKeKyKlb++HVwi6fqhVgG4xyw3GWPguzVAlBEPZDhDFiMQ9fnnn6\/+pLNIIbE1rEwkqjYWfHgh6xhznMJyxYcCcxyZ2SJZ9EZSs0tCdsYxmdeYrXW7RG\/MK6tz8Ma7H3Yxdqtx4gmiB3jQ5kDMICNkbNdEr0O+IlvQmDFjVJ5IhHyFZQQNNJzmmSpEb9yMxYaWHV3f6K9vzqRlR7oxv2jh7rHasEWOW6SkdFtaNmPfQNC33Xabmmu9evUK\/gkSDLxHIhG3vgHSA+YsImwiHR6J\/lgFQSTpBtdgJfTyyy+rEM+YM04UEr0TqIaoc\/PmzSoxAzYWjUQPAsHvYVVdccUV6k743UIagAZqJ4mDk0u+qOBxyKI3E72xTXZdN426vh2ij1buCeV\/D48LkJ6biR4EjaQm5jmkPwDGLFJWcwEhi2HNY4WDjcZw4ZKNm7HIgwp8unfvrqr2enJwvTrCRxKeNJDBYOAZE3prixkHqBCjCHMbuaiRehAHA50KNa1SCe7dq1yXrQoCI+JwF1aibpuznrDoAawV0eNrjnCzsaRlS2WirwzpQ0uPFIfHDtGD2LTnjW6HccM2VYgeAfXwIUXu0vr164e1DTTR46Jhw4apVYRfkoNjvwOrI0QXhaGGRDGIMvruu+8qyRGHphYsWKAShCNrGFwwkbEOcwgk72SGME30CHNtVbD3AjdPFBJ9VCbugYutiD5UomVYSLCy8HWF26FVrHHXuVcmwaK3Go5o\/fWBc7ikKKGeE4rojYSPtHoYSze8NJAIi4uL1capLunp6YKXO5R0g8TT8MSJZETAIkX\/sCJt3bp1xLfDKDfCXRbSZM+ePX1h0UfsvIjymwfRR8LVTl3RXKOlG3hYWZXvv\/9e7RVijNwwZ43t9LRFrze9AL72s9WdwwcAoIdzuTLukGMZnOxSddEi5V65D1ZBJdwrrdqPjWwQCvRL88lYu33WEg+00EhxeC6++OJgwhZdP+Q0EPtnn31WLiKn8fmwhqBRI6k3lsEIIfzTTz+54qXZuXOnDBkyRO3\/6AKihXwAjw\/IB3Y3Y\/X9+GhAgoB1CunBTvHrgSk7fcc1ySb6SPIYN2PtjqTFdaEs+ngRPR6pl8ExNjOm26stWSI1Ro2SHUOHyp62bWOqy3gziL60tFR5eVSW6I31gbTfeOMNFYcEP1YFMgR8yY877jhbfYE2O2nSJPUyu3kZbOxMKPdKfFQR6x5Zy0K5V+J+kDysfniGGbXl1atXqzMO+AlVUp3o9YEprdHbmlhxuIibsXEA0U4V0Ug3iHqH9HHwrzVb+sZnGZfB0DqPOuooO01x7BpY9JkjR8pfw4fHnehB9rm5uXEh+lCkD6ytgq8VFBSo4+zQUnF0XRe4GOK\/8TEwF7iCYokOnd5ty2BzW7GiBGmjH\/rAFHzrEXAOm4o4BQvCB1ngQA0MFHzQ8JFA7lf98cVKBh84fCSsTnimMtHDwNDhvzEGOFiGpEFOavN6rEn0jtFe+YqtiB4vBo6OczPWeiDCed3Ec\/gqo71Y7VcAAA54SURBVOubN2ON7YFG7xbpJhxOyEWLGCz42HXs2FFlC4N8dfTRRysvEBD4tddeq7R0JACHpAgZyKj36\/qxKdu7d2\/Lx9k9GRvPcWVd9mP+82RsjLMlGvdK6LtwuQRReMq9MkaMrG5PFNGbnx+tv745vr5XiN6hYQsr3ZgTj0DywcqJxRkENL6RNHrGuokRfyui1wemoBUDZH1gCstmHKAIF8TLdV43MWLkNqI3tgekD2kHZG4l8eB6+IXDowKWPom+4oiGmrOwIiEPebHANx0lHntHieh\/JKLnZmwMo4AgUJjgcFfDcXK4s0H71DomXPCgeWKyI9E23CqxRI4UqdHNiUdigKvCrcmy6K36YPeQFjYzt23b5nqNPp5jFakuvxkn8GZCbCBsPiczsUgk3O3+nURvF6kEXhdK70zg4xP2KLcRvdnSB\/GHiq8P694Lm7EJG0hDmOJE+5E71UcSvVPIlq\/XM370TsBBi94JVGOr06jrb926VUlxbve6ia3H0d3tuvhM0TW\/wtV+I3q3BuIj0fs0TLHxjXKzRR+OJ3S8GBL9AZT8tgr1G9HT6ybGL78Tt7tO73RRCAQn8I62TrfqndH2I57Xu27Oxtg5vxG9W+csLXpa9DG+qs7d7taXxrkeR66ZRB8Zo2Re4dY5S6IvLFRJSvCT9EKLvtwQuNUn2TxPKpMz1lgHgqUhAmu3bt2CWaes5iKJPulvadgGkOhdOD6u0ztJ9OVmiVv1TvNUjjZnrPF+fCRwqhthnZEIW6cXjET0rjFOYnyv\/SbduNU4oUXvpiQBJPpytOHW\/JvGRlYmZ6zxfsSiRxwcpMgcPHiwbaLv3LmzTJw4MUaaTf7tfiN6txonJHpq9Ml\/2y1a4NZlsLG5lckZq++HZDN69Gg5+eSTVZhjxGyya9HTj96d09atc5ZET6J35xtjiD3uZvfKyuaMhWSD2P7IEdCqVSuVUYlEn+XauWi3YSR6u0gl8DrXbWxRuik3+m59aYyNrGzOWEg2SJEHuWblypUkeoZAcJT5aNEXFqqY4LAaIxV8GJAoAnF2IhVEu0PuS2ipdqIKgtR2ffeddJ08WWTUqLAZpqKtG0lCsEkE69NOWyrTT2Bop0RTt1tPGcZK9JBskOAeGbQQv\/6XX36JmuiBN+I5IUyEVcEmMf7eoUOHsEODFHiI+KpTElpdjGtwLeoL91zcb\/fZeEcgWyHKbKQsb3brdKo\/6Fc4jIAPMpDhXXPbKpREX1io5rUdAgTB2r0W10UTOhbX5u7dK91E5KeCAtmUlhb25Yy2bifbHW3ddq\/XeLvhpYlXzlijZKNJozJEDwwjkS2IJ57XoS7UGem5+jonnu1EnfHqj8YbbXTDnDUSSEoTPazoW265xY4hymuSiABcCe2uGJxqZrxyxoaqx9jmxo0bh82Kxjnr1AjHt143zFkSfXzHlLWlMALR5owNtTJABqobbrhB+vfvL2eccYZKkZeI9HgpPGwp1\/WUtuhTbrTZYUcQiDZnLNILGks00o0jHWClvkeARO\/7IWYHnUYg2pyxZmudRO\/0CLF+Ej3nABEgAkTA5wiQ6H0+wOweESACRIBEzzlABIgAEfA5AiR6nw8wu0cEiAARINFzDhABIkAEfI4Aid7nA8zuEQE\/IIDzBzhvgBAkGRkZfuhSQvuQkkS\/Z88eee+991SMkfz8fFm6dKnKMIVJ5PXy559\/yqeffqr698ADD0i9evWCXYo1E5LXsfFy+708Z3Go7Nlnn5X69evL2rVr1Tt3ySWXSFaWdbRKfRANJ4F1ufLKK+XWW2919WEyvGOLFi2SRx55RMUyshMXKxHzMuWIHgPx9NNPy5dffiljxoyRunXrqkBEyPLzxBNPSMOGDROBuyPP+Pnnn+Wdd95R\/cHL9MILL5Qj+lgyITnSYFZqCwEvz1kQNoyoTp06qcBlf\/31l9x+++0qPPOgQYMsSXvGjBkydepUOfXUUxVGVatWlUMPPbTcfLYFXgIv2r59u7z55psyZcoUFdDQTTkDUo7oly9fLrAMLrroIrniiivUNFi1apWKHnjmmWeqyWg+uZjAuRLzo3QsFYS+NRJ9rJmQYm4YK6g0Al6eszA6Ro4cqeL3tGvXTmHw4YcfypAhQ2TcuHHSsWPHCrjs2LFD7rjjDvV+6nsqDV4SbtQ5Ckj0SQBfPxIWLxIxw4I\/\/vjj1a9LSkrUxNu2bZs8+eSTUqtWrSS2MLZHWxF9LJmQYmsR744VAa\/OWciIQ4cOlV9\/\/VWef\/55tcpEwcoThhXS7sGwMpe5c+dKv379lPWPMN+XXnqpdO3a1dWSjbEPJPpYZ3wc7r\/\/\/vvllVdeUUSv43TrCTlnzhxlZbRo0SIOT0pOFVZEX9lMSMnpBZ9qRMCrc3bz5s1q9YwgbUYDCpuqIPKjjz5aRo0aVWFz9bfffpN58+YpCQQSzu7du1X2Lfx4IdgbiT7J7y+0znvvvVclRcBSEincdMHLhA1M\/L5169ZJbmnlH29F9JXNhFT5lvDOeCDg5TlbVFQkV111lUqugkQrevNVfwCQAwK\/z8nJsYQKHwUkJoER9uijj0ZMjhIPzGOtg0QfK4Ix3u\/ll8Zu10n0dpHyxnVenrPxIHqMEvYoEML5uOOOUxu5brfqSfQueLdCLYOhBWLzZ+bMmWGTPrig+RGbEK10gyX1ww8\/7CoPgYidTLELvDpnraQb7fxw+OGHl7P0rYYVrqV33XWXbN261db1yZ4eJPpkj0DZjj+0vlTbjEUeTWxqYRl8\/vnnq5HQ1uLnn3+uPnBNmzZ1wQixCWYE4KXixTmrDaglS5ZEtRkbagbgY7d371657bbbXO8VR6J3wTscylUNuR6xadSjR4+Ucq+EGxt8mevUqaMsJp44dMEEDdEEL8\/ZUO6VMCzgW2\/lXmmGAPN08ODB0qtXL+ndu7c7B8nQKhK9C4ZIHz6ZNWuWytauD0xhowen91q2bOmCVla+CVbSDWqMlAmp8k\/lnU4i4OU5qw9MQV8fMGBA8MAUkrXAuMjMzJSJEyfKF198ofR3bNjCYeLggw9WK9Ds7GyZPn26Ou09bNgwqVGjhpNQx6VuEn1cYIy9ErhTTpo0SdasWSPNmjVTO\/rXXHON50m+uLhYvRB4gVDw4nTr1k19zFAiZUKKHVnW4BQCXp6z69atU0YUvG82bNigVo6XXXaZOq8C\/f2pp55SRsjYsWOVr\/3jjz8uIEskSj\/mmGPUO9q3b9+w3jlO4R5NvegLJFKEHsFZAIR5gNcRQqvgZG8yS8qdjE0m2Hw2ESACRCAZCJDok4E6n0kEiAARSCACJPoEgs1HEQEiQASSgQCJPhmo85lEgAgQgQQiQKJPINh8FBEgAkQgGQiQ6JOBOp9JBIgAEUggAiT6BILNRxEBIkAEkoEAiT4ZqPOZRIAIEIEEIkCiTyDYfBQRIAJEIBkIkOiTgTqfSQSIABFIIAIk+gSCzUcRASJABJKBAIk+GajzmUSACFRA4JdfflG5ZG+99VZPRKn00hCS6L00WmwrEfAxAghkhuBnL7zwggqAxhI\/BEj08cOSNREBIlBJBHR4bUS2ZG6ESoIY5jYSffwxZY1EgAhEiQBlmygBi\/JyEn2UgLntciRgRpq1t99+Wzp37ixDhgyRI444Qnbt2iVvvfWWShF4zz33yLHHHuu2prM9RCCIQDjZZvv27fLbb78JYvIjz2zNmjWJXJQIkOijBMyNl+tUcw0aNJAxY8YEE43g5VmxYoUUFha6Ps+mG3Flm2JHYOPGjfLJJ5\/IBx98IKtXr1aGB1JXjhgxQiZPniwTJkyQdu3aKQPFLNuA4JG8Hkk7kOpz6dKlKmcs6mjVqlXsjUuhGkj0PhhsZLaBVT9+\/Hh57rnnpEuXLgLNExmmLrjgAunYsaMPeskueBkB5IlFtiWkA0R6QeR+RTapnJwcZYyYvW2wUkUi+7y8vGCKwVdffVVloXriiSekQ4cOXoYj4W0n0ScccmceuGTJEunXr58ieWxmLVu2TB555BH1AdCpBJ15MmslApERWLVqlSJzpLtEqkDIjLqYZRsYLsjn\/M477ygPnDZt2kR+AK8IiwCJ3icTBPolLPgZM2aol2PevHmCpe91113nkx6yG15GoKSkRMkzmzZtUjlia9eurboTyttGb8wef\/zxSuKBpMMSGwIk+tjwc9XdSEgMq\/60005TG1dXX3210j9ZiECyEcA+Uv\/+\/aW0tLScxh7K2+bDDz+U66+\/Xkk3V1xxBfeX4jB4JPo4gOiWKnbs2CGDBw9Wm1ynnnqqjB49WmrUqOGW5rEdKYoA5iXkGkiIDz74oJJlevXqpax5bNKaD0lBysGmK\/R8EL0uMF6qV6+uNmdZokOARB8dXq6\/GpteAwYMUC\/JpZde6vr2soH+RGDfvn2yePFitdkKMofO3rp1a7Uhi41UOAhAvoFHjtnbZsGCBWru9ujRIyjdQOPHZixWBdxzin7OkOijx8zVd2ApfMstt8gDDzxAFzRXj5S\/G4dN15tvvlkWLVokQ4cOVXIiCuYlNlkRz6Z9+\/bKbdIc2wabse+9955yFW7SpImax3AdPv\/886VWrVr+Bs6h3pHoHQI2WdViMxZW0vDhw7mJlaxB4HOJgMsQING7bEBiaQ60UGxgnX322crNkoUIEAEiAARI9B6fBzh8cuONN6oj4s2aNVObrw899BCPiXt8XNl8IhBPBEj08UQzCXXBc2HkyJEqrk337t2VRQ9dk4UIEAEioBEg0XMuEAEiQAR8jgCJ3ucDzO4RASJABEj0nANEgAgQAZ8jQKL3+QCze0SACBABEj3nABEgAkTA5wiQ6H0+wOweESACRIBEzzlABIgAEfA5AiR6nw8wu0cEiAARINFzDhABIkAEfI4Aid7nA8zuEQEiQARI9JwDRIAIEAGfI0Ci9\/kAs3tEgAgQARI95wARIAJEwOcI\/D9gsD2DDl1b0gAAAABJRU5ErkJggg==","height":302,"width":378}}
%---
