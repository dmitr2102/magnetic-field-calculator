close all
% Выбираем точку наблюдения
obs_points = [0.05 0.05, 0.05];

% Создание катушки
coilX = Coil(0.09/2, 0.0002, 0.01, 45, 48, 100, [0 0 0], [1 0 0 90], 0.1);
coilX = coilX.generate();
[cBx1, cBx2, cBx3] = coilX.calculateField(obs_points);
cBx = [cBx1, cBx2, cBx3]

coilY = Coil(0.09/2, 0.0002, 0.01, 45, 48, 100, [0 0 0], [0 1 0 90], 0.1);
coilY = coilY.generate();
[cBy1, cBy2, cBy3] = coilY.calculateField(obs_points);
cBy = [cBy1, cBy2, cBy3]

coilZ = Coil(0.09/2, 0.0002, 0.01, 45, 48, 100, [0 0 0], [0 0 0 0], 0.1);
coilZ = coilZ.generate();
[cBz1, cBz2, cBz3] = coilZ.calculateField(obs_points);
cBz = [cBz1, cBz2, cBz3]

cB = cBx + cBy + cBz

figure
hold on
plot3(coilX.points(:,1),...
      coilX.points(:,2),...
      coilX.points(:,3), 'r', 'LineWidth', 0.001);
plot3(coilY.points(:,1),...
      coilY.points(:,2),...
      coilY.points(:,3), 'g', 'LineWidth', 0.001);
plot3(coilZ.points(:,1),...
      coilZ.points(:,2),...
      coilZ.points(:,3), 'b', 'LineWidth', 0.001);

valid = ~isnan(cB);
quiver3(obs_points(1), obs_points(2), obs_points(3), cB(1), cB(2), cB(3), 100, 'Color', 'b');

grid on
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);  % стандартный 3D-вид
xlim([-0.1, 0.1]);
ylim([-0.1, 0.1]);
zlim([-0.1, 0.1]);