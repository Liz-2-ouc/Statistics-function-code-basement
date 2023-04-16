function x = inv_terse(a,k)
    [m,n]=size(a);
    m_index=[1:m];%代表行
    n_index=[1:n];%代表列
    x(k,k)=1/a(k,k);
    x(k,n_index(n_index~=k))=a(k,n_index(n_index~=k))/a(k,k);
    x(m_index(m_index~=k),k)=-a(m_index(m_index~=k),k)/a(k,k);
    x(m_index(m_index~=k),n_index(n_index~=k))= ...,
        a(m_index(m_index~=k),n_index(n_index~=k))-a(m_index(m_index~=k),k)*a(k,n_index(n_index~=k))/a(k,k);
end