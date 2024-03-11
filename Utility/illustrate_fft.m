img=rgb2gray(imread('sherlock.jpg'));
A=fft(img,[],2);
R=20;

B=zeros(size(A));
B(:,1:R)=A(:,1:R);
B1=ifft(B,[],2);

E=zeros(size(A));
E(:,R:end/2)=A(:,R:end/2);
E1=ifft(E,[],2);

D=ifft(A-B-E,[],2);
% 
% E(:,100:end/2)=A(:,100:end/2);
% F=ifft(E,[],2);
% figure;subplot(1,3,1),imshow(A);
% subplot(1,3,2),imshow(B);
% subplot(1,3,3),imshow(A-B);

figure;subplot(2,2,1),imshow(img,[]);
subplot(2,2,3),imshow(B1,[]);
subplot(2,2,2),imshow(D,[]);
subplot(2,2,4),imshow(E1,[]);

