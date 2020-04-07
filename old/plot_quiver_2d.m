q1 = quiver(xvtg(:,:,1), yvtg(:,:,1), fields_interp(:,:,1,1), fields_interp(:,:,1,2), 0);
hold on;
q2 = quiver(xvtg(:,:,1), yvtg(:,:,1), fields_valid(:,:,1,1), fields_valid(:,:,1,2), 0);
hold off;

scale = 5;
qU1 = get(q1, 'UData');
qV1 = get(q1, 'VData');
set(q1, 'UData', scale*qU1, 'VData', scale*qV1);
qU2 = get(q2, 'UData');
qV2 = get(q2, 'VData');
set(q2, 'UData', scale*qU2, 'VData', scale*qV2);
