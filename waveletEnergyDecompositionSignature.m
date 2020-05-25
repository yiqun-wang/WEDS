function descriptor = waveletEnergyDecompositionSignature(shape_dir, shape_name, numEigs, weds_scale, ref, shape_ref, if_output, output_dir, collection_type)
%% 
addpath(genpath('utils'));

shape_fullname = [shape_dir, '/', shape_name];
if shape_name(end-2:end) == 'off'
    shape=read_shape(shape_fullname);
elseif shape_name(end-2:end) == 'obj'
    shape=read_shape_obj(shape_fullname);
elseif shape_name(end-2:end) == 'ply'
    shape=read_shape_ply(shape_fullname);
end

V = shape.VERT'; F = shape.TRIV';
[~, shape.n] = size(V);
options.symmetrize = 1;
options.normalize = 0;
type = 'conformal';
K = numEigs; 

[L,A] = compute_mesh_laplacian_plusA_half(shape.VERT',shape.TRIV',type,options);

try
    [V,D] = eigs(L, A, K+1, -1);
catch
    % In case of trouble make the laplacian definite
    [V,D] = eigs(L , A+ 1e-20*speye(shape.n), K+1, -1);
end
V=V(:,2:end); D=D(2:end,2:end);


if weds_scale <= 96
    sample = 3 + 2;
else
    sample = ceil(weds_scale/32) + 2;
end
if sample <= 32
    k_scale = floor(linspace(32,1,sample));
    k_scale = k_scale(2:end-1);
else
    k_scale = floor(linspace(32,1,sample-2));
end

Nscales=31;
if ref
    [L2,A2] = compute_mesh_laplacian_plusA_half(shape_ref.VERT',shape_ref.TRIV',type,options);
%     opts=struct('tol',5e-3,'p',10,'disp',0);
    lmax=max(eigs(L2,A2, K+1, -1));
else    
    lmax=max(diag(D));
end
lmax=lmax*1.01;
fprintf('Measuring largest eigenvalue, lmax = %g\n',lmax);
fprintf('Designing transform in spectral domain\n'); 
designtype='mexican_hat';
[g,t]=sgwt_filter_design(lmax,Nscales,'designtype',designtype);

Wnk=zeros(size(L,1),3, Nscales+1);
Tnk=zeros(size(L,1),size(L,1)); 
Win=zeros(length(k_scale),size(L,1),size(L,1));
clk=zeros(K,Nscales+1);
des=zeros(size(L,1),Nscales+1,3);

for k=1:numel(g)
  clk(:,k)=g{k}(diag(D));
  Wnk(:,:,k)=V*diag(clk(:,k))*V'*A*shape.VERT;
  Tnk=V*diag(clk(:,k))*V'*A;
  index = find(k_scale==k,1);  
  if ~isempty(index)
      if collection_type == 1
          WinT = abs(Tnk);
          Win(index,:,:) = ( (WinT-repmat(min(WinT),size(WinT,1),1))./repmat(max(WinT)-min(WinT),size(WinT,1),1) );
      elseif collection_type == 2
          WinT = Tnk.^2;
          Win(index,:,:) = (WinT./repmat(sqrt(sum(WinT.^2,1)),size(WinT,1),1));
      end
  end
end

for i=1:3
    wl=V'*A*squeeze(Wnk(:,i,:)).*clk;
    wl=diag(sum(wl,2));
    des(:,:,i)=(clk'*wl*(D.^2)*V'*A)'.*squeeze(Wnk(:,i,:));
end

if collection_type == 1
    descriptor=sum(des,3);
elseif collection_type == 2
    descriptor=sqrt(abs(sum(des,3)));
end

desc = [];
for i=1:length(k_scale)
    T = squeeze(Win(i,:,:));
    desc = [desc, (T'*descriptor)];  
end
descriptor = desc;
Fn = size(descriptor,2);
descriptor = descriptor(:,floor(linspace(1,Fn,weds_scale)));

if if_output
    Winn.V = V;
    Winn.clk = clk;
    Winn.A = diag(full(A));
    Winn.D = diag(D);
    save([output_dir, '/', shape_name(1:end-3), 'mat'],'-struct','Winn');
    dlmwrite([output_dir, '/', shape_name(1:end-3), 'txt'], descriptor,'delimiter',' ', 'precision','%4.6e');
end


end
