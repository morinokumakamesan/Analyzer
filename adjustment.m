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

Task_Num=xlsread('Task_frequency.csv','A1:A26');


%���f�[�^�������C�s��ɒǉ�
for i = 0:25
    if i == 0
        n = 0.1;
        m = i + 1;
    else
        n = i * 0.2;
        m = i + 1;
    end
    fa = fopen(['.\' num2str(Task_Num(m)) 'Hz.txt'],'r');    % �t�@�C���I�[�v���̕����̈����ݒ�͏o�͂���csv�t�@�C�����Ƃōl�����邱��
    dataA = fscanf(fa,'%f',[2,30000]);
    dataA = dataA.';
    fclose(fa);
    Data = vertcat(Data,dataA);

end

csvwrite('data.csv',Data,0,0);
