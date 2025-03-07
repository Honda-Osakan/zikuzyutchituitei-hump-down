# zikuzyutchituitei-hump-down
hump down

# zikuzyutchituitei-hump-down
clear; clc;

%% ===============================
%% 1. シミュレーション設定
%% ===============================
m1 = 1900;  m2 = 50;  m3 = 50;   % 車体、車軸の質量
k1 = 33e3;  k2 = 33e3;          % 前後サスペンション剛性 [N/m]
k3 = 260e3; k4 = 260e3;         % タイヤ剛性 [N/m]
c1 = 300;   c2 = 300;           % サスペンション減衰 [Ns/m]
l1 = 1.34;  l2 = 1.46;          % 重心から前・後車軸までの距離 [m]
J  = 3000;                      % 車体慣性モーメント [kg*m^2]

% 走行条件
v = 20;              % 車両速度 [m/s]
t_end = 10;          % シミュレーション時間 [s]
dt = 0.001;          
t = 0:dt:t_end;

% 状態変数の初期化
x1 = zeros(size(t));   dx1    = zeros(size(t));
theta = zeros(size(t)); dtheta = zeros(size(t));
x2 = zeros(size(t));   dx2    = zeros(size(t));
x3 = zeros(size(t));   dx3    = zeros(size(t));

Ff = zeros(size(t));  
Fr = zeros(size(t));
Fy1_array = zeros(size(t));  % 前車軸タイヤ力格納
Fy2_array = zeros(size(t));  % 後車軸タイヤ力格納

%% ===============================
%% 2. 台形型ハンプの設定
%% ===============================
% 台形型ハンプのパラメータ
%  ramp_up_length = 10m, top_length = 15m, ramp_down_length = 10m, hump_height = 0.2m
%  よってハンプは s=0(m) から s=35(m) まで続く
hump_start     = 0;     
ramp_up_length = 10;    
top_length     = 15;    
ramp_down_length = 10;  
hump_height    = 0.2;   

% ハンプの路面高さ関数
road_hump = @(s) ...
    (s < hump_start).*0 + ...
    (s >= hump_start & s < hump_start + ramp_up_length) .* ...
        ((s - hump_start)/ramp_up_length * hump_height) + ...
    (s >= hump_start + ramp_up_length & ...
     s < hump_start + ramp_up_length + top_length) .* ...
        hump_height + ...
    (s >= hump_start + ramp_up_length + top_length & ...
     s < hump_start + ramp_up_length + top_length + ramp_down_length) .* ...
        (hump_height - ...
         ((s - (hump_start + ramp_up_length + top_length)) ...
           / ramp_down_length * hump_height)) + ...
    (s >= hump_start + ramp_up_length + top_length + ramp_down_length).*0;

%% ===============================
%% 3. シミュレーションループ
%% ===============================
for i = 1:length(t)-1
    
    % 車両の現在進行距離 s [m]
    s = v * t(i);
    
    % 路面高さ (前後で同一と仮定)
    road_height = road_hump(s);
    y1 = road_height;
    y2 = road_height;
    
    % サスペンション力(前: Ff, 後: Fr)
    Ff(i) = -k1 * ((x1(i) - l1*theta(i)) - x2(i)) ...
            - c1 * ((dx1(i) - l1*dtheta(i)) - dx2(i));
    Fr(i) = -k2 * ((x1(i) + l2*theta(i)) - x3(i)) ...
            - c2 * ((dx1(i) + l2*dtheta(i)) - dx3(i));
    
    % タイヤ力 (前: Fy1, 後: Fy2)
    Fy1 = -k3 * (x2(i) - y1);
    Fy2 = -k4 * (x3(i) - y2);
    Fy1_array(i) = Fy1;
    Fy2_array(i) = Fy2;
    
    % 車体 (m1, J) の運動方程式
    ddx1    = 2*(Ff(i) + Fr(i)) / m1;
    ddtheta = 2*(-Ff(i)*l1 + Fr(i)*l2) / J;
    
    % 前軸 (m2), 後軸 (m3) の運動方程式
    ddx2 = ( -k3*(x2(i)-y1) - Ff(i) ) / m2;
    ddx3 = ( -k4*(x3(i)-y2) - Fr(i) ) / m3;
    
    % オイラー法で状態更新
    dx1(i+1)   = dx1(i)    + ddx1*dt;
    x1(i+1)    = x1(i)     + dx1(i)*dt;
    
    dtheta(i+1)= dtheta(i) + ddtheta*dt;
    theta(i+1) = theta(i)  + dtheta(i)*dt;
    
    dx2(i+1)   = dx2(i)    + ddx2*dt;
    x2(i+1)    = x2(i)     + dx2(i)*dt;
    
    dx3(i+1)   = dx3(i)    + ddx3*dt;
    x3(i+1)    = x3(i)     + dx3(i)*dt;
end

% 最終ステップの力を補完
Ff(end) = -k1*((x1(end)-l1*theta(end)) - x2(end)) ...
          - c1*((dx1(end)-l1*dtheta(end)) - dx2(end));
Fr(end) = -k2*((x1(end)+l2*theta(end)) - x3(end)) ...
          - c2*((dx1(end)+l2*dtheta(end)) - dx3(end));

final_s = v * t(end);
final_road = road_hump(final_s);
Fy1_array(end) = -k3*( x2(end) - final_road );
Fy2_array(end) = -k4*( x3(end) - final_road );

%% ===============================
%% 4. タイヤ力の可視化 (全区間)
%% ===============================
figure('Color','w');
plot(t, Fy1_array, 'b', 'LineWidth',1.5); hold on;
plot(t, Fy2_array, 'r', 'LineWidth',1.5);
xlabel('時間 [s]', 'FontSize', 16);
ylabel('垂直方向力 [N]', 'FontSize', 16);
legend({'前車軸タイヤ力','後車軸タイヤ力'}, 'FontSize', 16, 'Location','best');
grid on; set(gca,'FontSize',16);
title('');

%% ===============================
%% 5. ハンプ通過「後」の区間データを抽出
%% ===============================
% 全走行距離 s = v * t
s_all = v * t;

% ハンプは s = 0 から s = 35[m] で終了
% ここでは「ハンプ乗り終え」(s >= 35[m]) から、さらに一定区間 ( 35～50[m]) を計測期間とする
post_start = 35;   % ハンプ終了地点
post_end   = 50;   % その後 15m (約0.75秒) 走行分を測定対象とする

% インデックス抽出
post_idx = find(s_all >= post_start & s_all <= post_end);

% 抽出データ
t_post   = t(post_idx);           % 時間ベクトル (乗り終え後区間)
Fy1_post = Fy1_array(post_idx);   % 前車軸タイヤ力 (乗り終え後区間)
Fy2_post = Fy2_array(post_idx);   % 後車軸タイヤ力 (乗り終え後区間)

% 路面高さ
u_post   = arrayfun(@(ss) road_hump(ss), s_all(post_idx));

%% ===============================
%% 6. ARXモデルによる推定 (前車軸タイヤ力)
%% ===============================
% 前車軸タイヤ力 Fy1 に対する ARX モデル
data_Fy1_post = iddata(Fy1_post', u_post', dt);

% 例として [na=3, nb=3, nk=1] の ARX モデルを構築
na = 5; nb = 0; nk = 1;
model_arx_Fy1 = arx(data_Fy1_post, [na nb nk]);

% ARXモデルによる出力予測と評価
[simOut_Fy1, fitVal_Fy1, ~] = compare(data_Fy1_post, model_arx_Fy1);
y_pred_Fy1 = simOut_Fy1.OutputData;

% 誤差評価
err_Fy1   = Fy1_post' - y_pred_Fy1;
RMSE_Fy1  = sqrt(mean(err_Fy1.^2));
maxErr_Fy1= max(abs(err_Fy1));

%% ===============================
%% 7. ARXモデルによる推定 (後車軸タイヤ力)
%% ===============================
% 後車軸タイヤ力 Fy2 に対する ARX モデル
% 図示はしないが、評価は出力する
data_Fy2_post = iddata(Fy2_post', u_post', dt);

model_arx_Fy2 = arx(data_Fy2_post, [na nb nk]);
[simOut_Fy2, fitVal_Fy2, ~] = compare(data_Fy2_post, model_arx_Fy2);
y_pred_Fy2 = simOut_Fy2.OutputData;

err_Fy2   = Fy2_post' - y_pred_Fy2;
RMSE_Fy2  = sqrt(mean(err_Fy2.^2));
maxErr_Fy2= max(abs(err_Fy2));

%% ===============================
%% 8. 結果表示
%% ===============================
fprintf('\n=== ハンプ通過「後」の区間における ARX モデル適用結果 ===\n\n');

disp('--- 前車軸タイヤ力 Fy1 ---');
disp(model_arx_Fy1);
fprintf('Fit率         = %.2f %%\n', fitVal_Fy1);
fprintf('RMSE          = %.4f [N]\n', RMSE_Fy1);
fprintf('最大ピーク誤差 = %.4f [N]\n\n', maxErr_Fy1);

disp('--- 後車軸タイヤ力 Fy2 (図示なし) ---');
disp(model_arx_Fy2);
fprintf('Fit率         = %.2f %%\n', fitVal_Fy2);
fprintf('RMSE          = %.4f [N]\n', RMSE_Fy2);
fprintf('最大ピーク誤差 = %.4f [N]\n', maxErr_Fy2);

%% ===============================
%% 9. 前車軸タイヤ力予測の比較プロット
%% ===============================
figure('Color','w');
plot(t_post, Fy1_post, 'b-', 'LineWidth',1.5); hold on;
plot(t_post, y_pred_Fy1, 'm--', 'LineWidth',2.5);
xlabel('時間 [s]', 'FontSize', 16);
ylabel('垂直方向力 [N]', 'FontSize', 16);
legend({'モデルによる実測出力','モデルによる出力予測'}, ...
       'FontSize', 16, 'Location','best');
grid on; set(gca,'FontSize',16);
title('');
hump down

