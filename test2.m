m = 2;
n = 2;
m1 = 201;
n1 = 201;
tpixel = max(m1,n1);
NAfily = 69.7077;
NAfilx = 69.7077;
NApixel = 2*max(round(NAfily),round(NAfilx));
x = linspace(-tpixel/NApixel,tpixel/NApixel,tpixel);
[X,Y] = meshgrid(x,x);
[theta,r] = cart2pol(X,Y);
idx = r<=1;
z = zeros(size(X));

r =r(idx);
theta = theta(idx);

m_abs = abs(m);
rpowers = [];
length_r = length(r);

for j = 1:length(n)
    rpowers = [rpowers m_abs(j):2:n(j)];
end
rpowers = unique(rpowers);

if rpowers(1)==0
    rpowern = arrayfun(@(p) r.^p,rpowers(2:end),'UniformOutput',false);
    rpowern = cat(2,rpowern{:});
    rpowern = [ones(length_r,1) rpowern];
else
    rpowern = arrayfun(@(p)r.^p,rpowers,'UniformOutput',false);
    rpowern = cat(2,rpowern{:});
end

z = zeros(length_r,length(n));
for j = 1:length(n)
    s = 0:(n(j)-m_abs(j))/2;
    pows = n(j):-2:m_abs(j);

    for k = length(s):-1:1
        p = (1-2*mod(s(k),2))* ...
                   prod(2:(n(j)-s(k)))/              ...
                   prod(2:s(k))/                     ...
                   prod(2:((n(j)-m_abs(j))/2-s(k)))/ ...
                   prod(2:((n(j)+m_abs(j))/2-s(k)));
        idx = (pows(k)==rpowers);
        z(:,j) = z(:,j) + p*rpowern(:,idx);
    end
end
idx_pos = m>0;
idx_neg = m<0;

if any(idx_pos)
    z(:,idx_pos) = z(:,idx_pos).*cos(theta*m_abs(idx_pos)');
end
if any(idx_neg)
    z(:,idx_neg) = z(:,idx_neg).*sin(theta*m_abs(idx_neg)');
end
zOut = zeros(size(X));
zOut = z;