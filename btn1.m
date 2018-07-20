function btn1

global text
global input_num
global input_num2
global input_num3
global input_num4
global input_num5
global input_num6
global time
global input_freq
global input_freq2
global input_freq3
global input_freq4
global input_freq5
global input_freq6
global input_loop
global loop_check
global minForce

Task = zeros(0,0,0);
Precision = zeros(0,0);
CV = zeros(0,0);
vpi_AE = zeros(0,0);
pvi_AE = zeros(0,0);
VPI_Precision = zeros(0,0);
PVI_Precision = zeros(0,0);
minForce = 1;

%loopTaskの解析
if loop_check.Value == 1
    Task_Num=xlsread('Task_frequency.csv','A1:A26');
    
    %26tasksというディレクトリがなければ作成する
    if exist('26tasks', 'file') == 0
        mkdir 26tasks;
    end
    
    %input_loop(標準は26回)まわしてdata.csvを各周波数ごとのcsvにわける
    for i = 1:str2num(input_loop.String)
        display(['Analyzing...' num2str(i) '/' input_loop.String]);
        time = (str2num(input_num.String) * str2num(input_loop.String)) * 1000;
        AllData=xlsread('data.csv',['A1:B' num2str(time) '']);
        ActualData=xlsread('data.csv',['B1:B' num2str(time) '']);
        
        Split_AllData = AllData((i-1)*30000+1:30000*i,:);
        Split_ActualData = ActualData((i-1)*30000+1:30000*i,:);
        
        period = Task_Num(i) * str2num(input_num.String);
        
        %初期化
        Adjustment_Data = zeros(0,0);
        Adjustment_All = zeros(0,0);
        vpi_AE = zeros(0,0);
        pvi_AE = zeros(0,0);
        
        %周期ごとにhorzcatする Adjustment_Dataの最初の列にはターゲットが
        for j = 1:period
            Adjustment_All = horzcat(Adjustment_All, Split_AllData((j-1)*fix(1000/Task_Num(i))+1:fix(1000/Task_Num(i))*j,:));
            Adjustment_Data = horzcat(Adjustment_Data, Split_ActualData((j-1)*fix(1000/Task_Num(i))+1:fix(1000/Task_Num(i))*j,:));
        end

        %【PVIVPI】
        [m,n] = size(Adjustment_Data);
        
        for p = 1:period
            for q = 1:round(m/2)
                vpi_AE(q+round(m/2)*(p-1),1) = abs(Adjustment_All(q,2*p-1) - Adjustment_All(q,2*p));
                pvi_AE(q+round(m/2)*(p-1),1) = abs(Adjustment_All((q-1)+round(m/2),2*p-1) - Adjustment_All((q-1)+round(m/2),2*p));
            end
        end
        
        vpi_Target = Adjustment_Data(1:round(m/2),1);
        pvi_Target = Adjustment_Data(round(m/2):m,1);
        
        VPI_Target = sum(vpi_Target);
        PVI_Target = sum(pvi_Target);
        VPI_Target = VPI_Target * (n-1);
        PVI_Target = PVI_Target * (n-1);
        
        %【/PVIVPI】
        
        min = zeros(0,0);
        max = zeros(0,0);
        %総データ数回数まわす Split_AllData(k,3)に各絶対誤差
        for k = 1:str2num(input_num.String)*1000
            Split_AllData(k,3) = abs(Split_AllData(k,1) - Split_AllData(k,2));
            
            %【CV】
            if Split_AllData(k,1) == 4
                max = horzcat(max,Split_AllData(k,2));
            elseif Split_AllData(k,1) == 1
                min = horzcat(min,Split_AllData(k,2));
            end
            %【/CV】
        end
        
        %【Precision】
        %絶対誤差
        All_Sum = sum(Split_AllData);
        All_AE = All_Sum(1,3);
        
        VPI_AE = sum(vpi_AE);
        PVI_AE = sum(pvi_AE);
        
        disp(All_Sum(1,1));
        disp(VPI_Target + PVI_Target);
        disp(m*(n-1));
        disp(minForce * round(m/2) * (n-1));
        disp(All_Sum(1,1) - minForce * str2num(input_num.String)*1000);
        disp((VPI_Target - minForce * round(m/2) * (n-1)) + (PVI_Target - minForce * round(m/2) * (n-1)));
        disp((VPI_Target - minForce * round(m/2) * (n-1)));
        disp(All_AE);
        disp(VPI_AE + PVI_AE);
        disp(VPI_AE);
        disp(PVI_AE);
        disp(size(vpi_AE));
        disp(size(pvi_AE));
        
        %Precisionの算出
        precision = ((All_Sum(1,1) - minForce * str2num(input_num.String)*1000) - All_AE) / (All_Sum(1,1) - minForce * str2num(input_num.String)*1000) * 100;
        result_precision = horzcat(Task_Num(i),precision);
        Precision = vertcat(Precision, result_precision);
        
        precision_vpi = ((VPI_Target - minForce * round(m/2) * (n-1)) - VPI_AE) / (VPI_Target - minForce * round(m/2) * (n-1)) * 100;
        result_precision_vpi = horzcat(Task_Num(i),precision_vpi);
        VPI_Precision = vertcat(VPI_Precision, result_precision_vpi);
        
        precision_pvi = ((PVI_Target - minForce * round(m/2) * (n-1)) - PVI_AE) / (PVI_Target - minForce * round(m/2) * (n-1)) * 100;
        PVI_Precision = vertcat(PVI_Precision, precision_pvi);
        %【/Precision】

        %【CV】
        cv_min = std(min) / mean(min) * 100;
        cv_max = std(max) / mean(max) * 100;
        cv = horzcat(cv_max,cv_min);
        result_cv = horzcat(Task_Num(i), cv);
        CV = vertcat(CV, result_cv);
        %【/CV】
        
        %26Taskを周期ごとに分けた結果を保存
        %csvwrite(['26tasks/' num2str(Task_Num(i)) 'Hz.csv'],Adjustment_Data,0,0);
    end
    VPIPVI = horzcat(VPI_Precision,PVI_Precision);
    xlswrite('Precision.xlsx',Precision,1,'A1');
    xlswrite('VPIPVI.xlsx',VPIPVI,1,'A1');
    xlswrite('CV.xlsx',CV,1,'A1');

%6Taskの解析
else
    %--------------------------------------------工事中-------------------------------------------------
    time = (str2num(input_num.String) + str2num(input_num2.String) + str2num(input_num3.String) + str2num(input_num4.String) + str2num(input_num5.String) + str2num(input_num6.String)) * 1000;
    AllData=xlsread('data.csv',['A1:B' num2str(time) '']);
    ActualData=xlsread('data.csv',['B1:B' num2str(time) '']);
    display(['A1:A' num2str(time) '']);
    %6回のTaskを各シートに分ける
    for i =1:6
        Task(:,:,i) = AllData((i-1)*30000+1:30000*i,:);
        xlswrite('6Task_data.xlsx',Task(:,:,i),i,'A1');
    end
end

%各タスクのデータを周期ごとにまとめる


%xlswrite('data_write.xlsx',Period(:,:,1),1,'A1');
display('Analysis completed!');
end