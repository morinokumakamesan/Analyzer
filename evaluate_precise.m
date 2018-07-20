Data = zeros(0,0);
AData = zeros(0,0);
num = zeros(0,0);
Pre = zeros(0,0);
SDataMin = zeros(0,0);
SDataMax = zeros(0,0);
min = zeros(0,0);
max = zeros(0,0);
MinSum = zeros(0,0);
MaxSum  = zeros(0,0);
Rise = zeros(0,0);
Fall = zeros(0,0);

%生データを順次，行列に追加
for i = 0:25
    if i == 0
        n = 0.1;
    else
        n = i * 0.2;
    end
    fa = fopen(['.\' num2str(n) 'Hz.txt'],'r');    % ファイルオープンの部分の引数設定は出力したcsvファイルごとで考慮すること
    dataA = fscanf(fa,'%f',[2,30000]);
    dataA = dataA.';
    fclose(fa);
    Data = horzcat(Data,dataA);
    for j = 1:30000
        Data(j,(i+1)*3) = abs(Data(j,(i)*3+1) - Data(j,(i)*3+2));
    end
end

%csvwrite('sample.csv',Data,0,0);

%ブレの総和，精度の導出
for i = 0:25
    if i == 0
        nl = 0.1;
    else
        nl = i * 0.2;
    end
    num = horzcat(num,nl);
    sum = 0;
    for j = 1:30000
        sum = sum + Data(j,(i+1)*3);
    end
    AData = horzcat(AData,sum);
    pre = ((45000 - sum) / 45000) * 100;
    Pre = horzcat(Pre,pre);
end

num = num.';
Pre = Pre.';
AData = AData.';

Out = horzcat(num,Pre);
Out = horzcat(Out,AData);

%csvwrite('ALL.csv',Out,0,0);


for i = 0:25
    if i == 0
        nl = 0.1;
    else
        nl = i * 0.2;
    end
    min = zeros(0,0);
    max = zeros(0,0);
    sum_min = 0;
    sum_max = 0;
    f = 0;
    g = 0;
    for j = 1:30000
        if Data(j,i*3+1) == 1
            min = horzcat(min,Data(j,i*3+2));
            f = f + 1;
        elseif Data(j,i*3+1) == 4
            max = horzcat(max,Data(j,i*3+2));
            g = g + 1;
        end
    end
    minsum = 0;
    maxsum = 0;
    k = 0;
    l = 0;
    for k = 1:f
        minsum = minsum + min(1,k);
    end
    minsum = std(min) / (minsum / f) *100;
    for l = 1:g
        maxsum = maxsum + max(1,l);
    end
    maxsum = std(max) / (maxsum / g) * 100;
    MinSum = horzcat(MinSum,minsum);
    MaxSum = horzcat(MaxSum,maxsum);
end

MinSum = MinSum.';
MaxSum = MaxSum.';

CV = horzcat(num,MinSum);
CV = horzcat(CV,MaxSum);

%csvwrite('CV.csv',CV,0,0);


%VPI PVI

for i = 0:25
    if i == 0
        nl = 0.1;
    else
        nl = i * 0.2;
    end
    %num = horzcat(num,nl);
    rise = zeros(0,0);
    fall = zeros(0,0);
    f = 0;
    g = 0;
    rd = 0;
    fd = 0;
    RD = zeros(0,0);
    FD = zeros(0,0);
    for j = 1:30000-1
        if Data(j,i*3+1) < Data(j+1,i*3+1)
            %rd = Data(j,i*3+2) - Data(j,i*3+1); %M-T
            %rise = horzcat(rise,rd);
            rise = horzcat(rise,Data(j,i*3+3));
            f = f + 1;
        elseif Data(j,i*3+1) > Data(j+1,i*3+1)
            %fd = Data(j,i*3+2) - Data(j,i*3+1); %M-T
            %fall = horzcat(fall,fd);
            fall = horzcat(fall,Data(j,i*3+3));
            g = g + 1;
        end
    end
    R = 0;
    F = 0;
    m = 0;
    n = 0;
    for m = 1:f
        R = R + rise(1,m);
    end
    %R = R / f;
    R = ((22500 - R) / 22500) * 100;
    for n = 1:g
        F = F + fall(1,n);
    end
    %F = F / g;
    F= ((22500 - F) / 22500) * 100;
    Rise = horzcat(Rise,R);
    Fall = horzcat(Fall,F);
    
    if i == 2
        disp(rise);
        disp(fall);
        csvwrite('rise.csv',rise,0,0);
        csvwrite('fall.csv',fall,0,0);
    end
end

Rise = Rise.';
Fall = Fall.';

VPIPVI = horzcat(num,Rise);
VPIPVI = horzcat(VPIPVI,Fall);

csvwrite('VPIPVI.csv',VPIPVI,0,0);


Time = 30;
SampleN = Time * 1000;
CN = zeros(0,0);
ind = zeros(0,0);
fr = 0;
datan = 0;
for i = 0:25;  %周波数[Hz]
    if i == 0
        fr = 0.1;
    else 
        %datan = i * 0.2;
        fr = i * 0.2;
    end
count = 0;
division = 0;
fl = 0;
DD = zeros(0,0);
Ind = zeros(0,0);
Ind2 = zeros(0,0);
PeakN = fr * Time;
division = SampleN / PeakN;
m = 0;
n = 0;
o = 0;
p = 0;
for j = 1:30000
    count = j - fl * division;
    %DD = horzcat(DD,Data(j,i*3+2));
    if count <= division
    else
        %disp (numel(Ind));
        %disp (numel(Ind2));
        [m,n] = size(Ind);
        [o,p] = size(Ind2);
        if n < p
            Ind2(:,p) = [];
        end
        Ind = vertcat(Ind,Ind2);
        count = 1;
        fl = fl + 1;
        Ind2 = zeros(0,0);
    end
    CN = horzcat(CN,count);
    if fl == 0
        Ind = horzcat(Ind,Data(j,i*3+2));
    else
        Ind2 = horzcat(Ind2,Data(j,i*3+2));
    end
end
[m,n] = size(Ind);
[o,p] = size(Ind2);
if n == p
else
    Ind2(:,p) = [];
end
Ind = vertcat(Ind,Ind2);
%Ind = vertcat(Ind,mean(Ind));
%Ind = vertcat(Ind,std(Ind));

Ind = Ind.';
%CNOUT = horzcat(CN,Ind);
%csvwrite(['CN_' num2str(fr) '.csv'],Ind,0,0);
end