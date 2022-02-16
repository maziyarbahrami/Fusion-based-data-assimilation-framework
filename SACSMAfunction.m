%FUNCTON SAC - SMA
%m:input vector contains 21 columns:
%parameters (13 ones),state variables (6 ones), precipitation and pet
%parameters definition
%p1:uztwm.upper zone tension water capacity
%p2:uzfwm:upper zone free water capacity
%p3:lztwm:lower zone tension water capacity
%p4:lzfpm:lower zone free primary capacity
%p5:lzfsm:lower zone free supplementary water capacity
%p6:adimp:additional impervious area
%p7:uzk:upper zone depletion parameter
%p8:lzpk:lower zone prtimary depletion parameter
%p9:lzsk:lower zone seconadary depletion parameter
%p10:zperc:maximum percolation rate
%p11:rexp:percolation equation exponent
%p12:pctim:impervious area of watershed
%p13:pfree:free water percolation from upper to lower zone
%%
%state variables definitions
%s1:utwc:upper zone tension water content
%s2:uzfwc:upper zone free water content
%s3:lztwc:lower zone tension water content
%s4:lzfpc:lower zone free primary water cotent
%s5:lzfsc:lower zone free secondary water content
%s6:adimc:additional impervious area water content
%20th column:precipitation
%21th:pet
function [y1,y2,y3,s1,s2,s3,s4,s5,s6]=sacsma(m)
tz=0.00001;
        parea=1-m(6)-m(12);
        et1=m(21)*(m(14)/m(1));
        %         et1(:,:);
        red=m(21)-et1;
        m(14)=m(14)-et1;
        %         m(k,14)
%         red(k,t)

        et2=0;
        if m(14)<=0
            et1=et1+m(14);
            m(14)=0;
            red=m(21)-et1;
            
            if m(15)<red
                et2=m(15);
                m(15)=0;
                red=red-et2;
                if m(14)<tz
                    m(14)=0;
                end
                if m(15)<tz
                    m(15)=0;
                end
            else
                et2=red;
                m(15)=m(15)-et2;
                red=0;
            end
        else
            if (m(14)/m(1))<(m(15)/m(2))
                uzrat=(m(14)+m(15))/(m(1)+m(2));
                m(14)=m(1)*uzrat;
                m(15)=m(2)*uzrat;
            end
            if m(14)<tz
                m(14)=0;
            end
            if m(15)<tz
                m(15)=0;
            end
        end
        
%        m(k,15) 
        et3=red*m(16)/(m(1)+m(3));
        m(16)=m(16)-et3;
        if m(16)<0
            et3=et3+ m(16);
            m(16)=0;
        end
        et5=et1+(red+et2)*(m(19)-et1-m(14))/(m(1)+m(3));
        m(19)=m(19)-et5;
        if m(19)<0
            et5=et5+m(19);
            m(19)=0;
        end
        et5=et5*m(6);

        twx=m(20)+m(14)-m(1);
%         twx(k,t)
        if twx<0
            m(14)=m(14)+m(20);
            twx=0;
        else
            m(14)=m(1);
        end
        m(19)=m(19)+m(20)-twx;
        roimp=m(20)*m(12);
        sbf=0;
        ssur=0;
        sif=0;
        sperc=0;
        sdro=0;
        ninc=floor(1+0.2*(m(15)+twx));
%         ninc(k,t)
%         m(k,15)
        dinc=1/ninc;
        pinc=twx/ninc;
%         pinc(k,t)
        duz=1-(1-m(7))^dinc;
%         duz(k,t)
        dlzp=1-(1-m(8))^dinc;
        dlzs=1-(1-m(9))^dinc;
        ninc=1;
        for n=1:ninc
%             addro=0;
            bf_p=0;
            bf_s=0;
%             percm=0;
%             perc=0;
%             defr=0;
%             adsur=0;
%             ratio=0;
%             check=0;
%             del=0;
%             perct=0;
%             percf=0;
%             hpl=0;
%             ratlp=0;
%             ratls=0;
%             fracp=0;
%             percp=0;
%             percs=0;
%             fracp=0;
%             sur=0;
            %%
            adsur=0;
            ratio=(m(19)-m(14))/m(3);
%             ratio(k,t,n)
            if ratio<0
                ratio=0;
            end
            addro=pinc*ratio^2;
            bf_p=m(17)* dlzp;
%             bf_p(k,t)
            m(17)=m(17)-bf_p;
            if m(17)<=0.0001
                bf_p=bf_p+m(17);
                m(17)=0;
            end
            sbf=sbf+bf_p;
%             sbf(k,t)
            bf_s=m(18)*dlzs;
%             bf_s(k,t)
            m(18)=m(18)-bf_s;
            if m(18)<0.0001
                bf_s=bf_s+m(18);
                m(18)=0;
            end
            sbf=sbf+bf_s;
            if (pinc+m(15))<=0.01
                m(15)=m(15)+pinc;
            else
                percm=m(4)*dlzp+m(5)*dlzs;
%                 percm(k,t)
                perc=percm*m(15)/m(2);
%                 perc(k,t)
                defr=1-(m(16)+m(17)+m(18))/(m(3)+m(4)+m(5));
%                 defr(k,t)
                if defr<0
                    defr=0;
                end
                perc=perc*(1+m(10)*(defr^m(11)));
                if perc>=m(15)
                    perc=m(15);
                end
                m(15)=m(15)-perc;
                check=m(16)+m(17)+m(18)+perc-m(3)-m(4)-m(5);
                if check>0
                    perc=perc-check;
                    m(15)=m(15)+check;
                end
                sperc=sperc+perc;
%                 sperc(k,t)
                del=m(15)*duz;
%                 del(k,t)
                sif=sif+del;
                m(15)=m(15)-del;
                perct=perc*(1-m(13));
%                 perct(k,t)
                if (perct+m(16))<=m(3)
                    m(16)=m(16)+perct;
                    percf=0;
                else
                    percf=m(16)+perct-m(3);
                    m(16)=m(3);
                    
                end
                percf=percf+(perc*m(13));
                if percf~=0
                    hpl=m(4)/(m(4)+m(5));
                    ratlp=(m(17)/m(4));
                    ratls=(m(18)/m(5));
                    fracp=hpl*2*(1-ratlp)/(2-ratlp-ratls);
                    if fracp>1
                        fracp=1;
                        
                    end
                    percp=percf*fracp;
                    percs=percf-percp;
                    m(18)= m(18)+percs;
                    if  m(18)>m(5)
                        percs=percs-m(18)+m(5);
                        m(18)=m(5);
                        
                    end
                    m(17)=m(17)+percf-percs;
                    if m(17)>=m(4)
                        excess=m(17)-m(4);
                        m(16)=m(16)+excess;
                        m(17)=m(4);
                    end
                end
                if pinc~=0
                    if ((pinc+m(15))<=m(2))
                        m(15)=m(15)+pinc;
                    else
                        sur=pinc+m(15)-m(2);
%                         sur(k,t)
                        m(15)=m(2);
                        ssur=ssur+(sur*parea);
                        adsur=sur*(1-addro/pinc);
                        ssur=ssur+adsur*m(6);
                    end
                    
                    
                end
                
            end
%             m(k,19)
            m(19)=m(19)+pinc-addro-adsur;
%             pinc(k,t)
%             m(k,19)
            if m(19)>(m(1)+m(3))
                addro=addro+ m(19)-(m(1)+m(3));
                m(19)=m(1)+m(3);
            end
            sdro=sdro+(addro*m(6));
            if m(19)<tz
                m(19)=0;
            end
            
            
            
        end
        eused=et1+et2+et3;
%         eused(k,t)
        sif=sif*parea;
%         sif(k,t)
        tbf=sbf*parea;
        base=tbf;
        surf=roimp+sdro+ssur+sif;
        eused=eused*parea;
        tet=eused+et5;
%         tet(k,t)
        if m(19)<m(14)
            m(19)=m(14);
        end
        tot_outflow=surf+base;
        if tot_outflow<0
            tot_outflow=0;
            surf=0;
            base=0;
        end
        simflow= tot_outflow;
        surf_tot=surf;
        base_tot=base;
        s25=sqrt(0.15*surf_tot/3);
%         y=surf_tot+s25*randn;
y1=surf;
y2=base;
y3=y1+y2;
%OR TOT-OUTFLOW??????????????????????????????
        s1=m(14);
        s2=m(15);
        s3=m(16);
        s4=m(17);
        s5=m(18);
        s6=m(19);
end
