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
global check_1
global check_2
global check_3
global check_4
global check_5
global check_6

Task = zeros(0,0,0);
Precision = zeros(0,0);
CV = zeros(0,0);
vpi_AE = zeros(0,0);
pvi_AE = zeros(0,0);
VPI_Precision = zeros(0,0);
PVI_Precision = zeros(0,0);
minForce = 1;
Precision2 = zeros(0,0);
stack = 0;

%loopTaskの解析
if loop_check.Value == 1
    Task_Num=xlsread('Task_frequency.csv','A1:A26');
    
    %26tasksというディレクトリがなければ作成する
    if exist('26tasks', 'file') == 0
        mkdir 26tasks;
    end
    if exist('26-SinglePeriods', 'file') == 0
        mkdir 26-SinglePeriods;
    end
    
    %input_loop(標準は26回)まわしてdata.csvを各周波数ごとのcsvにわける
    for i = 1:str2num(input_loop.String)
        display(['Analyzing...' num2str(i) '/' input_loop.String]);
        time = (str2num(input_num.String) * str2num(input_loop.String)) * 1000;
        AllData=xlsread('data.csv',['A1:B' num2str(time) '']);
        ActualData=xlsread('data.csv',['B1:B' num2str(time) '']);
        TargetData=xlsread('data.csv',['A1:A' num2str(time) '']);
        
        %上から，各タスクごとに切り出したターゲットとセンサ値（AllData），センサ値のみ（ActualData），ターゲットのみ（TargetData）
        %Split_AllData = AllData((i-1)*30000+1:30000*i,:);
        Split_ActualData = ActualData((i-1)*30000+1:30000*i,:);
        Split_TargetData = TargetData((i-1)*30000+1:30000*i,:);
        
        %ローパスフィルタ（カットオフ10Hz）
        [b,a]=butter(4,.02);
        Split_ActualData = filtfilt(b,a,Split_ActualData);
        
        Split_AllData = horzcat(Split_TargetData,Split_ActualData);        
        
        period = Task_Num(i) * str2num(input_num.String);
        
        %初期化
        Adjustment_Data = zeros(0,0);
        Adjustment_All = zeros(0,0);
        vpi_AE = zeros(0,0);
        pvi_AE = zeros(0,0);
        vpi_Target = zeros(0,0);
        pvi_Target = zeros(0,0);
        Adjustment_Target = zeros(0,0);
        adjustment = zeros(0,0);
        add = 1;
        
        %1周期目のVPIの"1"の個数と最後の周期PVIの"1"の個数を数える（50サンプルの中から抽出）
        B = size(find(Split_TargetData(1:50,1) == 1));
        C = size(find(Split_TargetData(str2num(input_num.String)*1000-49:str2num(input_num.String)*1000,1) == 1));
        disp(B);
        disp(C);
        
        %周期ごとにhorzcatする Adjustment_Dataの最初の列にはターゲットが
        for j = 1:period
            
            %1周期のデータ数に変動がない場合は，そのまま各周期ごとに切る
            if rem(1000,Task_Num(i)) == 0
                Adjustment_All = horzcat(Adjustment_All, Split_AllData((j-1)*fix(1000/Task_Num(i))+1:fix(1000/Task_Num(i))*j,:));
                Adjustment_Data = horzcat(Adjustment_Data, Split_ActualData((j-1)*fix(1000/Task_Num(i))+1:fix(1000/Task_Num(i))*j,:));
                Adjustment_Target = horzcat(Adjustment_Target, Split_TargetData((j-1)*fix(1000/Task_Num(i))+1:fix(1000/Task_Num(i))*j,:));
            else
                %最後の周期以外
                if j ~= period
                    %以下の条件に当てはまった場合は"0"を追加，当てはまらなかった場合は"次"のデータを追加
                    %切り分けた周期の最後が"1"かどうか
                    O = Split_TargetData(j*fix(1000/Task_Num(i))+add,1) == 1.0;
                    %本来あるべき"1"の数（B）に収まっているか
                    P = size(find(Split_TargetData((j-1)*fix(1000/Task_Num(i))+add:(j-1)*fix(1000/Task_Num(i))+add+50,1)==1)) >= B(1,1)-1;
                    Q = size(find(Split_TargetData((j-1)*fix(1000/Task_Num(i))+add:(j-1)*fix(1000/Task_Num(i))+add+50,1)==1)) <= B(1,1)+1;
                    R = size(find(Split_TargetData(j*fix(1000/Task_Num(i))+add:j*fix(1000/Task_Num(i))+add+50,1)==1)) >= C(1,1)-1;
                    S = size(find(Split_TargetData(j*fix(1000/Task_Num(i))+add:j*fix(1000/Task_Num(i))+add+50,1)==1)) <= C(1,1)+1;
                    
                    %いろんな条件にあてはまったとき
                    if O(1,1) && (P(1,1) && Q(1,1)) && (R(1,1) && S(1,1))
                        adjustment = vertcat(Split_AllData((j-1)*fix(1000/Task_Num(i))+add:fix(1000/Task_Num(i))*j+(add-1),:),zeros(1,2));
                        adjustment_target = vertcat(Split_AllData((j-1)*fix(1000/Task_Num(i))+add:fix(1000/Task_Num(i))*j+(add-1),1),0);
                        adjustment_data = vertcat(Split_AllData((j-1)*fix(1000/Task_Num(i))+add:fix(1000/Task_Num(i))*j+(add-1),2),0);
                    else
                        adjustment = vertcat(Split_AllData((j-1)*fix(1000/Task_Num(i))+add:fix(1000/Task_Num(i))*j+(add-1),:),Split_AllData(j*fix(1000/Task_Num(i))+1,:));
                        adjustment_target = vertcat(Split_TargetData((j-1)*fix(1000/Task_Num(i))+add:fix(1000/Task_Num(i))*j+(add-1),1),Split_TargetData(j*fix(1000/Task_Num(i))+add,1));
                        adjustment_data = vertcat(Split_AllData((j-1)*fix(1000/Task_Num(i))+add:fix(1000/Task_Num(i))*j+(add-1),2),Split_AllData(j*fix(1000/Task_Num(i))+add,2));
                        add = add + 1;
                    end
                %最後の周期
                else
                    adjustment = Split_AllData((j-1)*fix(1000/Task_Num(i))+add:str2num(input_num.String)*1000,:);
                    adjustment_target = Split_AllData((j-1)*fix(1000/Task_Num(i))+add:str2num(input_num.String)*1000,1);
                    adjustment_data = Split_AllData((j-1)*fix(1000/Task_Num(i))+add:str2num(input_num.String)*1000,2);
                    
                    T = size(Split_AllData((j-1)*fix(1000/Task_Num(i))+add:str2num(input_num.String)*1000,:));
                    if T(1,1) == fix(1000/Task_Num(i))
                        adjustment = vertcat(Split_AllData((j-1)*fix(1000/Task_Num(i))+add:str2num(input_num.String)*1000,:),zeros(1,2));
                        adjustment_target = vertcat(Split_AllData((j-1)*fix(1000/Task_Num(i))+add:str2num(input_num.String)*1000,1),0);
                        adjustment_data = vertcat(Split_AllData((j-1)*fix(1000/Task_Num(i))+add:str2num(input_num.String)*1000,2),0);
                    else
                        adjustment = Split_AllData((j-1)*fix(1000/Task_Num(i))+add:str2num(input_num.String)*1000,:);
                        adjustment_target = Split_AllData((j-1)*fix(1000/Task_Num(i))+add:str2num(input_num.String)*1000,1);
                        adjustment_data = Split_AllData((j-1)*fix(1000/Task_Num(i))+add:str2num(input_num.String)*1000,2);
                    end
                    
                end
                Adjustment_All = horzcat(Adjustment_All, adjustment);
                Adjustment_Data = horzcat(Adjustment_Data, adjustment_data);
                Adjustment_Target = horzcat(Adjustment_Target, adjustment_target);
            end
            
        end

        %--------------------PVIVPI--------------------
        [m,n] = size(Adjustment_Data);
        
        %苦渋の決断（1周期の最後のデータはあったりなかったりするので，含めてしまうとVPIに比べPVIの絶対誤差は低くなりがちとなってしまう）
        m=m-1;
        
        for p = 1:period
            for q = 1:round(m/2)
                vpi_AE(q+round(m/2)*(p-1),1) = abs(Adjustment_All(q,2*p-1) - Adjustment_All(q,2*p));
                pvi_AE(q+round(m/2)*(p-1),1) = abs(Adjustment_All((q-1)+round(m/2),2*p-1) - Adjustment_All((q-1)+round(m/2),2*p));
                vpi_Target(q+round(m/2)*(p-1),1) = Adjustment_All(q,2*p-1);
                pvi_Target(q+round(m/2)*(p-1),1) = Adjustment_All((q-1)+round(m/2),2*p-1);
            end
        end
        
        VPI_Target = sum(vpi_Target);
        PVI_Target = sum(pvi_Target);
        
        %--------------------/PVIVPI--------------------
        
        min = zeros(0,0);
        max = zeros(0,0);
        %総データ数回数まわす Split_AllData(k,3)に各絶対誤差
        for k = 1:str2num(input_num.String)*1000
            Split_AllData(k,3) = abs(Split_AllData(k,1) - abs(Split_AllData(k,2)));
            
            %--------------------CV--------------------
            if Split_AllData(k,1) == 4
                max = horzcat(max,Split_AllData(k,2));
            elseif Split_AllData(k,1) == 1
                min = horzcat(min,Split_AllData(k,2));
            end
            %--------------------/CV--------------------
        end
        
        %--------------------Precision--------------------
        %絶対誤差
        All_Sum = sum(Split_AllData);
        All_AE = All_Sum(1,3);
        
        VPI_AE = sum(vpi_AE);
        PVI_AE = sum(pvi_AE);
        
        %Precisionの算出
        precision = ((All_Sum(1,1) - minForce * str2num(input_num.String)*1000) - All_AE) / (All_Sum(1,1) - minForce * str2num(input_num.String)*1000) * 100;
        result_precision = horzcat(Task_Num(i),precision);
        Precision = vertcat(Precision, result_precision);
        
        precision_vpi = ((VPI_Target - minForce * round(m/2) * (n)) - VPI_AE) / (VPI_Target - minForce * round(m/2) * (n)) * 100;
        result_precision_vpi = horzcat(Task_Num(i),precision_vpi);
        VPI_Precision = vertcat(VPI_Precision, result_precision_vpi);
        
        precision_pvi = ((PVI_Target - minForce * round(m/2) * (n)) - PVI_AE) / (PVI_Target - minForce * round(m/2) * (n)) * 100;
        PVI_Precision = vertcat(PVI_Precision, precision_pvi);
        
        precision2 = (2*PVI_Target - 2 * minForce * round(m/2) * (n) - (PVI_AE+VPI_AE)) / (2*PVI_Target - 2 * minForce * round(m/2) * (n)) * 100;
        result_precision2 = horzcat(Task_Num(i),precision2);
        Precision2 = vertcat(Precision2, result_precision2);
        %--------------------/Precision--------------------

        %--------------------CV--------------------
        cv_min = std(min) / mean(min) * 100;
        cv_max = std(max) / mean(max) * 100;
        cv = horzcat(cv_max,cv_min);
        result_cv = horzcat(Task_Num(i), cv);
        CV = vertcat(CV, result_cv);
        %--------------------/CV--------------------
        
        %--------------------SinglePeriod--------------------
        %飛び出た分の個数をかぞえる
        rema = str2num(input_num.String)*1000 - fix(1000/Task_Num(i)) * n;
        
        if rem(1000,Task_Num(i)) == 0
            SinglePeriod = horzcat(mean(Adjustment_Target,2),mean(Adjustment_Data,2));
        else
            at_all = Adjustment_Target(fix(1000/Task_Num(i))+1,:);
            ad_all = Adjustment_Data(fix(1000/Task_Num(i))+1,:);

            Adjustment_Target_fix = Adjustment_Target(1:fix(1000/Task_Num(i)),:);
            Adjustment_Data_fix = Adjustment_Data(1:fix(1000/Task_Num(i)),:);
            SinglePeriod = horzcat(mean(Adjustment_Target_fix,2),mean(Adjustment_Data_fix,2));
            %飛び出た分の付け足し
            SinglePeriod = vertcat(SinglePeriod,horzcat(sum(at_all,2)/rema,sum(ad_all,2)/rema));
        end
        %--------------------/SinglePeriod--------------------
        
        %26Taskを周期ごとに分けた結果を保存
        csvwrite(['26tasks/' num2str(Task_Num(i)) 'Hz.csv'],Split_AllData,0,0);
        xlswrite(['26tasks/' num2str(Task_Num(i)) 'Hz.xlsx'],Adjustment_Target,1,'A1');
        xlswrite(['26-SinglePeriods/' num2str(Task_Num(i)) 'Hz.xlsx'],SinglePeriod,1,'A1');
    end
    VPIPVI = horzcat(VPI_Precision,PVI_Precision);
    xlswrite('Precision.xlsx',Precision,1,'A1');
    xlswrite('VPIPVI.xlsx',VPIPVI,1,'A1');
    xlswrite('CV.xlsx',CV,1,'A1');

    
%6Taskの解析
else
    %--------------------------------------------工事中-------------------------------------------------
    check_ar = [check_1.Value check_2.Value check_3.Value check_4.Value check_5.Value check_6.Value];
    time_ar = [str2num(input_num.String) str2num(input_num2.String) str2num(input_num3.String) str2num(input_num4.String) str2num(input_num5.String) str2num(input_num6.String)];
    
    valid_data = check_ar.*time_ar;

    AllData = xlsread('data.csv',['A1:B' num2str(sum(valid_data)*1000) '']);
    ActualData = xlsread('data.csv',['B1:B' num2str(sum(valid_data)*1000) '']);
    
    for k = 1:sum(valid_data)*1000
        AllData(k,3) = abs(AllData(k,1) - abs(AllData(k,2)));
    end
    
    disp(sum(check_ar));
    %6回のTaskを各シートに分ける
    for i =1:sum(check_ar)
        xlswrite('6Task_data.xlsx',AllData(stack+1:(valid_data(i)*1000)+stack,:),i,'A1');
        stack = stack + valid_data(i)*1000;
        disp(stack);
    end
end

%各タスクのデータを周期ごとにまとめる

%xlswrite('data_write.xlsx',Period(:,:,1),1,'A1');
display('Analysis completed!');
end