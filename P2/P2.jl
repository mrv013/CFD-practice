### A Pluto.jl notebook ###
# v0.19.2

using Markdown
using InteractiveUtils

# ╔═╡ 0e70748b-8b49-44a2-9318-2938e727421d
using LinearAlgebra

# ╔═╡ c2351b4f-89df-4a33-a073-92413245464c
using SparseArrays

# ╔═╡ 1e30de32-7004-4498-867c-7d8e028a4ee7
using Plots 

# ╔═╡ b0f3964c-9f15-4148-904f-815e0d78b017
"""<div style="direction: rtl;text-align:center"> <h1>  بنام خدا </h1></div>
<div style="direction: rtl;text-align:center">
<body>

<body></div>
"""|>HTML

# ╔═╡ 77b0c19c-0b5a-4b62-9b56-eb97b6434e2a
"""<div style="direction: rtl;text-align:center">
<body>
 <p><b>دانشگاه صنعتی امیر کبیر</b> <p>
 <p><b>تکلیف سری اول: مطالعه انتقال حرارت یک بعدی در فین</b> <p>
 <p><b>عنوان درس:</b>
CFD
 <p>
 <p> <b>دانشجو:</b>
محمد رضا واعظی  <p>
 <p><b>استادر راهتما:</b>
دکتر نادران 
 <p>
 <p><b>تاریخ تحویل:</b>

1400/12/26
<p>

<body></div>
"""|>HTML

# ╔═╡ 7ce121eb-ff77-4c54-8d8f-07d286d603f0
"""<div style="direction: rtl"> <h3>  مقدمه و پارامتر ها </h3><hr </div>
<style>

*{box-sizing:border-box;}

body,div,h2{
 font face: "Times New Roman"
 size :"5"
 direction: rtl;
 padding:0;
 margin:0;

}

.my-header{
 direction: rtl;
 padding:15px;
 #text-align:center;
}

.my-content{
 font face = "Comic sans MS"
 #font face = "Times New Roman"
 #font face: "Verdana"
 size :"5"
 direction: rtl;
 padding:14px;
 #text-align:center;
}

.my-footer{
 padding:15px;
 text-align:center;
}

</style>
<div class=my-content> <p>در این بخش پارامتر های مسئله را تعریف میکنیم همان طور که  از فرمول مشخص است ما 4 نوع ثابت در نظر گرفته ایم شامل : K (ضریب رسانایی حرارت)، hc (ضریب همرفتی )و c1,c2 که ضرایب ساده کننده در فرمول هستند </p> </div>

"""|>HTML

# ╔═╡ c529e753-f54a-4eb4-aa23-3ec975e27a84


# ╔═╡ 56eee02f-091a-4b5c-a850-08fb59f12aae
 GC.gc()

# ╔═╡ 9a1a0057-54a5-4c33-93ff-c3cdd70f1662
k=0.1

# ╔═╡ eaa783cb-f687-4075-b320-3f7050636e04
hc=0.01

# ╔═╡ 7e5e2c1a-6b51-43fd-83e4-78cd416615ff
c1=1

# ╔═╡ ac541cbb-c35c-490b-a0a6-4daef8cd3566
c3=1

# ╔═╡ aece47b6-9aa9-474d-ab82-da97eb53c2ad
"""

<head>
<style>

*{box-sizing:border-box;}

body,div,h2{
 font face: "Times New Roman"
 size :"5"
 direction: rtl;
 padding:0;
 margin:0;

}

.my-header{
 direction: rtl;
 padding:15px;
 #text-align:center;
}

.my-content{
 font face = "Comic sans MS"
 #font face = "Times New Roman"
 #font face: "Verdana"
 size :"5"
 direction: rtl;
 padding:14px;
 #text-align:center;
}

.my-footer{
 padding:15px;
 text-align:center;
}font face = "Times New Roman" size :"5
<head>
</style>
<div style="direction: rtl;padding:14px";
 #text-align:center;> <p>tau , h به ترتیب فاصله زمانی و مکانی مسئله هتند و در این بخش به تعریف آن ها اقدام میکنیم  </p> </div>

"""|>HTML

# ╔═╡ ee228b61-0538-49f0-ab30-74246458cf51
τ=0.01

# ╔═╡ 5be587ad-7703-4430-a9eb-279a24d20fc0
"""
<div div style="direction: rtl;padding:14px"> <p> برای پیدا کردن h باید ابتدا تعداد نقاط مورد نظر برای حل مسئله را مشخص کنیم سپس با تابع تعریف شده مقدار h را تعیین میکنیم  </p> </div>

"""|>HTML

# ╔═╡ 68b3dff3-1a99-4c31-ad3d-3b2455de69a6
n=5

# ╔═╡ b90a98a0-e085-4f36-b415-c89a0bff1624
dh(n)=1/(n-1)

# ╔═╡ 5a4bc4e6-1bad-4a4a-bd73-00d56055d895
h=dh(n)

# ╔═╡ 60e0b122-25da-42a5-80e9-bf42bb6f6455
"""
<div div style="direction: rtl;padding:14px"> <p>در پایان زمان پایانی برای حل مسئله را مشخص و مقادیر اولیه دما و دمای محیط را مشخص میکنیم </p> </div>

"""|>HTML

# ╔═╡ b5eeb804-2cc8-436c-9acc-1277eb664195
end_time=3

# ╔═╡ 1b92dbcf-92c3-42e7-b7cb-c60bd5af3642
T_inf=20

# ╔═╡ 59cd42b9-5792-429a-8d56-5e007d693150
T0=100

# ╔═╡ e4c74479-e565-4570-bd3a-fc4ec8c52b46
"""<div style="direction: rtl"> <h2>  حل عددی </h2><hr </div>


"""|>HTML

# ╔═╡ b4f62d03-a975-4ce4-9de1-f6fd12944982
"""
<div div style="direction: rtl;padding:14px"> <p> اضافه کردن کتاب خانه های مورد استفاده در کد را در این بخش انجام میدهیم    </p> </div>

"""|>HTML

# ╔═╡ 2201cf8f-8e48-47db-a5f0-f6cc2265ef8f
"""
<div div style="direction: rtl;padding:14px"> <p>برای این که در در هر مرحله متغیر جدید تعریف نکنیم و زمان اجرای کد را به حداقل برسانیم متعیر هایی که میخوایم را با ایت تابع تعریف اولیه میکنیم (اصافه کردن و اصلاح حافظه جز فرآیند های زمان بر است و بهتر است به یکباره انجام شود)</p> </div>

"""|>HTML

# ╔═╡ 3b7e65df-07de-4c5c-8e42-e7159aeeb752
function initilize_variables(n)
	return (
    Vector{Int64}(undef,3*(n-2)+3),
	Vector{Int64}(undef,3*(n-2)+3),
	Vector{Float64}(undef,3*(n-2)+3),
	zeros(Float64,n),
	zeros(Float64,n)
	)
end

# ╔═╡ 1f3db9fb-0156-4f7d-bd6c-14fd4f512a45
"""
<div div style="direction: rtl;padding:14px"> <p> در ادامه برای هر چهار روش FTCS، BTCS، CTCS و کرنک جدا گانه توابع و بررسی خطا وابسته به بازه های زمانی و مکانی را انجام داده و بررسی آن ها را تکمیل میکنیم .  </p> </div>

"""|>HTML

# ╔═╡ 866d2f41-fdc3-4e49-ba91-49d95999eeed
"""<div style="direction: rtl"> <h2> انجام محاسبات با روش FTCS</h2><hr </div>


"""|>HTML

# ╔═╡ 76164d36-37e4-4d6b-8d7b-b315c1c5c801
"""
<div div style="direction: rtl;padding:14px"> <p>  هر روش محاسبه 4 نوع تابع نیاز دارد اول باید ماتریس ظرایب را به روش OCC تعریف میکنیم سپس تابعی برای معرفی شرایط مرزی تعریف میکنیم و بعد به سراغ دو تابع بکب برای انجام محاسبات در دامنه مکان با یک قدم در زمان و دیگری برای گسترش محاسبات تا end _time  استفاده میکنیم    </p>
<p> تابع تعریف ضرایب برای نقات میانی :</p>
</div>

"""|>HTML

# ╔═╡ be2591de-b5f1-4e16-b7d9-56c864c362db
function internal_coeffs_FTCS!(I, J, V, C, n, h, τ, c1, k, c3, T_inf )
	for i= 2:(n-1)
		K=3*(i-2).+(1:3)# [3*n+1 ][3*o+1, 3*o+2, 3*o+3]
		I[K] .= i#[i, i, i]
		J[K] = [i-1, i, i+1]
		V[K] = [k*τ/c1/h^2, (1.0 -2*k*τ/c1/(h^2)-c3*τ/c1), k*τ/c1/(h^2)]
		C[i] = c3*τ*T_inf/c1
	end
	#return I, J, V, C
end

# ╔═╡ 9bd7a2c0-5ab2-4781-bd25-6cf42edc0e96
"""
<div div style="direction: rtl;padding:14px"> 
<p> تابع تعریف ضرایب در نقات مرزی :</p>
</div>
"""|>HTML

# ╔═╡ c5c42ec3-463a-4814-b9b3-eac861f1b7bc
function left_bc!(I, J, V, C, n)
        K= 3*(n-2) .+(1:1)
		#k=3*(i-2).+(1:3)
		I[K] .= 1
		J[K] = [1]
		V[K] = [1]
	    C[1] = 0
	   # b[1] = T0
	#return I, J, V, C
end

# ╔═╡ e9bfb747-0921-4cc7-ab6b-20ef741c3453
"""
<div div style="direction: rtl;padding:14px"> <p> دقت شود تابع شرایط مرزی در شمت راست باید حتما بعد از اجرا تابع ضرایب باشد چرا که بخشی از ضرایب مرزی با اصلاح ضرایب میانی تکمیل شده است </p>
<p> تابع تعریف ضرایب در نقات مرزی سمت راست :</p>
</div>
"""|>HTML

# ╔═╡ b31236ce-3469-4759-9697-26ffa86566fc
function right_bc!(I, J, V, C, n, h, hc, k, T_inf)
        K= 3*(n - 2)+1 .+(1:2)
	    p=3*(n-3)
		I[K] .= n
		J[K] = [n - 1, n]
	    V[K] = V[p+1:p+2]+ V[p+3]*[-2*hc*h/k, 1]
	    #V[K] = [k*τ/c1/h^2, (1.0 -2*k*τ/c1/(h^2)-c3*τ/c1)]+k*τ/c1/(h^2)*[-2*hc*h/k 1]
		#V[K] =(k/(k+hc*h))*V[p]
	    C[n] = (V[p+3])*2*T_inf*hc*h/k + C[n-1]
	return I, J, V, C
end

# ╔═╡ 7f1cd37b-27fb-45fb-961e-9bfd2e979b55
"""
<div div style="direction: rtl;padding:14px">  
<p> تعریف تابع برای محاسبه تابع Dense:</p>
</div>
"""|>HTML

# ╔═╡ aea48b90-9784-46ec-83b8-e6757b428184
	bc = (left_bc!,right_bc!)

# ╔═╡ fc63e35c-0b20-4a36-8760-543db4f0cac5
function generate_les_FTCS(n, hc, k, τ, c1, c3, T_inf,bc)
    h=dh(n)
	lbc,rbc=bc
	I, J, V, C= initilize_variables(n)
	internal_coeffs_FTCS!(I, J, V, C, n, h, τ, c1, k, c3, T_inf)
	lbc(I, J, V, C, n)
	rbc(I, J, V, C, n, h, hc, k, T_inf)
	return sparse(I, J, V) , C
end

# ╔═╡ c1cbc128-b47d-40b2-a969-ddfe919ca11c
"""
<div div style="direction: rtl;padding:14px">  
<p> تابع تعریف ماتریس شرایط اولیه(در زمان صفر):</p>
</div>
"""|>HTML

# ╔═╡ 9298374c-c9af-4cef-ab61-e53922636d75
time0=[T0;ones(n-1,1)*T_inf]

# ╔═╡ a68b1cca-7642-4593-8ba0-04b36931c20a
A_FTCS ,C_FTCS = generate_les_FTCS(n, hc, k, τ, c1, c3, T_inf,bc)

# ╔═╡ 70ea54df-7d89-4411-90a4-56230db9c868
det(A_FTCS)

# ╔═╡ 3a950fd4-b631-4868-bcb7-07089d2eb4f9
"""
<div div style="direction: rtl;padding:14px">  
<p>  تعریف تابع برای محاسبه در کشتره زمان (همان طور که مشخص است یکی از ورودی های این تابع end_time  است ):</p>
</div>
"""|>HTML

# ╔═╡ 0b754307-067b-4ab8-b857-7b340c82eb50
function generate_all_time_FTCS(n, hc, k, τ, c1, c3,T0,T_inf,end_time,bc)
	h_range=LinRange(0,1,n)
    time0=[T0;ones(n-1,1)*T_inf]
	A ,C = generate_les_FTCS(n, hc, k, τ, c1, c3, T_inf,bc)
	time_lapse=collect(τ:τ:end_time)
	cloum=length(time_lapse)
	all_time=zeros(Float64,n,cloum)
	for i in (1:cloum)
		T1=A*(time0)+C
		all_time[:,i]=T1
		time0=T1
	end
	return (all_time_FTCS=all_time , time_lapse_FTCS=time_lapse, h_range_FTCS=h_range)
end

# ╔═╡ e2649731-93ce-4783-942f-82b6954bf061
T1_FTCS=generate_all_time_FTCS(n, hc, k, τ, c1, c3,T0,T_inf,end_time,bc)

# ╔═╡ b2398175-5d89-4311-9651-0b494984a37a
plo1s=scatter(T1_FTCS.h_range_FTCS,T1_FTCS.all_time_FTCS[:,1],xlabel="ξ",ylabel="θ", lable="n=$n")

# ╔═╡ befb0951-9cc9-4f5c-ac6d-81ae08842b66
"""
<div style="direction: rtl"> <h3>   تحلیل خطا در روش FTCS</h3><hr </div>



<div div style="direction: rtl;padding:14px"> <p> برای تحیلی خطا به 3 توع عمل میکنیم
<b>نوع اول: </b>
با استفاده از توابع تعریف شده مسئله را برای چنداز n (تعداد نقاط در طول فین)با مقدار ثابت از tau انجام  و مقایسه مکنیم. 
<b>نوع دوم: </b>
خطای sucsecive  را برای نوع اول در نقاط مشخص محاسبه میکنیم.
<b>نوع سوم: </b>
  محاسبات را برای چند بازه زماتی(tau) و با تغییر   n  انجام میدهیم و باهم مقایه میکنیم.    
   </p>



<p> </p>
</div>
"""|>HTML

# ╔═╡ fdc80205-b52b-4606-a6f9-16288c674299
nRange1=[1,2,4,8].*(n-1).+1

# ╔═╡ 1878e5d0-ff15-4a6d-a2da-30f86bf9ab22
length(nRange1)

# ╔═╡ e218cf86-2e2b-412f-b6ed-e5d6ba340f39
"""
<div div style="direction: rtl;padding:14px">  
<p> محاسبات برای n موجود در مجموعه nRange:</p>
</div>
"""|>HTML

# ╔═╡ f68d5d7e-3f64-4b04-9aab-31448a59fd3c
res_FTCS= [ generate_all_time_FTCS(ni, hc, k, τ, c1, c3,T0,T_inf,end_time,bc) for ni in nRange1 ]

# ╔═╡ f069dfeb-0cc6-49df-86c3-5aa60f807f99
"""
<div div style="direction: rtl;padding:14px">  
<p> نمایش خطای نوع اول برای n موجود در مجموعه nRange:</p>
</div>
"""|>HTML

# ╔═╡ 583a8f2c-913b-4eda-ad90-f66b3a38d87b
begin
	
	nr1=length(nRange1)
    plotc=scatter(xlabel="h",ylabel="FTCS in time=0.5", title="slution for diffrent gride")
    for i in 1:nr1-1 
	  scatter!(plotc, res_FTCS[i].h_range_FTCS, res_FTCS[i].all_time_FTCS[:,100], label="n=$(nRange1[i])" )
    end
	plotc
end

# ╔═╡ 7ff114c1-ebcd-48b8-8e5e-ec0e2168dfed
"""
<div div style="direction: rtl;padding:14px">  
<p> همان طور که مشخص است برای n=33  روش FTCS دچار نا پایداری میشود که در نمودار زیر مشخص است  </p>
</div>
"""|>HTML

# ╔═╡ cb06c89c-b0fc-43d4-b145-5724bfcbfb63
plotx=scatter(res_FTCS[4].h_range_FTCS, res_FTCS[4].all_time_FTCS[:,100], label="n=$(nRange1[4])",xlabel="h",ylabel="FTCS in time=0.5", title="slution for diffrent gride" )

# ╔═╡ 3689c5df-404d-483d-91a1-c3df537deb62
"""
<div div style="direction: rtl;padding:14px">  
<p> در این مرحله برای محاسبه خطای successive ابتدا باید نقاط ثابتی را در طول فین انتخاب کنیم و شماره نقاط را در هر n پیدا کرده و سپس به محاسبه خطا و رسم نمودار پرداخت . تابع تعریف شده index  به این منظور است .  </p>
</div>
"""|>HTML

# ╔═╡ 9e47327b-2e81-4b48-9331-09a5e2aba870
index2(ξ,Δξ)=Int(ξ/Δξ)+1

# ╔═╡ c2a528ad-c689-4ce9-84dd-4fe532cb317d
index2.(res_FTCS[1].h_range_FTCS,res_FTCS[2].h_range_FTCS[2])

# ╔═╡ 0a9f70c5-9be2-44d3-8541-7f0111977ca6
begin
	err=zeros(nr1-1,n)
	h3= @. 1.0/(nRange1-1)
	plo3=scatter(axis= :log,xlabel="\$ h\$",ylabel="successive Error", title="slution convergence")
    for i in 1:( nr1-1 )
        i1=index2.(res_FTCS[1].h_range_FTCS,res_FTCS[i].h_range_FTCS[2])
	    i2=index2.(res_FTCS[1].h_range_FTCS,res_FTCS[i+1].h_range_FTCS[2])
		for j in 2:nRange1[1]
	    err[i,j]= abs.(res_FTCS[i+1].all_time_FTCS[i2[j],100]-res_FTCS[i].all_time_FTCS[i1[j],100])
		end
    end
	lable=["ξ=$ξ" for i in 1:1,ξ in res_FTCS[1].h_range_FTCS[2:end] ]
	scatter(plo3, h3[1:end-1], err[:,2:end], label=lable )
end


# ╔═╡ c3da6280-8d6e-4e26-94d2-e1cd3a3e292b
begin
	
	#nr1=length(nRange1)
    plott=scatter(xlabel="h",ylabel="FTCS times in end of fin ", title="slution for diffrent gride",legend= :bottomright)
    for i in 1:nr1-1 
	  scatter!(plott, res_FTCS[i].time_lapse_FTCS, res_FTCS[i].all_time_FTCS[end,:], label="n=$(nRange1[i])" )
    end
	plott
end

# ╔═╡ b7684934-c28f-4d12-bab5-7b521adfe9bb
"""<div style="direction: rtl"> <h2>   انجام محاسبات با روش BTCS</h2><hr </div>

<div div style="direction: rtl;padding:14px">  
<p> از انجایی که مطالب توضیح داده شده در روش محاسبه FTCS بسیاری با سایر روش ها مشابه هستند تنها به ذکر مواتفاوت ها میپردازیم .  </p>
</div>


"""|>HTML

# ╔═╡ d1357278-e475-4e85-8491-306847558306
function internal_coeffs_BTCS!(I, J, V, C, n, h, τ, c1, k, c3, T_inf)
	for i= 2:(n-1)
		K=3*(i-2).+(1:3)# [3*n+1 ][3*o+1, 3*o+2, 3*o+3]
		I[K] .= i#[i, i, i]
		J[K] = [i-1, i, i+1]
		V[K] = [-k*τ/c1/h^2, (1.0 +2*k*τ/c1/h^2+c3*τ/c1), k*τ/c1/h^2]
		C[i] = -c3*τ*T_inf/c1
	end
end

# ╔═╡ 8d192a41-7f8e-48ea-ba0c-7bd94872e1f8
function generate_les_BTCS(n, hc, k, τ, c1, c3, T_inf,bc)
    h=dh(n)
	lbc,rbc=bc
	I, J, V, C= initilize_variables(n)
	internal_coeffs_BTCS!(I, J, V, C, n, h, τ, c1, k, c3, T_inf)
	lbc(I, J, V, C, n)
	rbc(I, J, V, C, n, h, hc, k, T_inf)
	return sparse(I, J, V) , C
end

# ╔═╡ 55b27c1d-8ae1-49b6-986a-eab21ea2383d
A_BTCS ,C_BTCS = generate_les_BTCS(n, hc, k, τ, c1, c3, T_inf,bc)

# ╔═╡ 71390fa2-e25b-4039-89ef-9f915e0b68d3
det(A_BTCS)

# ╔═╡ 8af39678-beaf-4d30-8d21-7704c8601219
function generate_all_time_BTCS(n, hc, k, τ, c1, c3,T0,T_inf,end_time,bc)
    h_range=LinRange(0,1,n)
	time0=[T0;ones(n-1,1)*T_inf]
	A ,C = generate_les_BTCS(n, hc, k, τ, c1, c3,T_inf,bc)
	time_lapse=collect(τ:τ:end_time)
	cloum=length(time_lapse)
	all_time=zeros(Float64,n,cloum)
	for i in (1:cloum)
		T1=A\((time0)-C)
		all_time[:,i]=T1
		time0=T1
	end
	return (all_time=all_time , time_lapse=time_lapse, h_range=h_range)
end


# ╔═╡ 7880a698-2149-46bf-9010-b55bf32d22fb
generate_all_time_BTCS(n, hc, k, τ, c1, c3, T0, T_inf, end_time,bc)

# ╔═╡ 280e11f6-4edb-449d-a873-91ed7161dabb
res_BTCS= [ generate_all_time_BTCS(ni, hc, k, τ, c1, c3,T0,T_inf,end_time,bc) for ni in nRange1 ]

# ╔═╡ 78808f18-5bc5-4ee4-825b-eadcf3c3ca45
all_time_BTCS=generate_all_time_BTCS(n, hc, k, τ, c1, c3,T0,T_inf,end_time,bc)

# ╔═╡ 69fdbe90-7ce4-4a09-8b45-018d425d741b
"""
<div style="direction: rtl"> <h3>   تحلیل خطا در روش BTCS</h3><hr </div>



<div div style="direction: rtl;padding:14px"> <p> از آنجا که روش BTCS روش ضمنی است خطا برای مانند روش قبلی (FTCS) در هیچ مقدار از n دچار ناپیداری نمیشود و این موضوع را در نمودار ها میتوان دید .


   </p>



<p> </p>
</div>
"""|>HTML

# ╔═╡ 928a6042-3c76-47c3-8c76-7caf7397b981
begin
	
	#nr1=length(nRange1)
    plotBTCS_n=scatter(xlabel="h",ylabel="BTCS in time=0.5", title="slution for diffrent gride")
    for i in 1:nr1 
	  scatter!(plotBTCS_n, res_BTCS[i].h_range, res_BTCS[i].all_time[:,100], label="n=$(nRange1[i])" )
    end
	plotBTCS_n
end

# ╔═╡ 5aa48d02-15a0-4747-8dfe-c4cc76c86fbc
"""
<div div style="direction: rtl;padding:14px">  
<p> نمایش خطای نوع دوم برای n موجود در مجموعه nRange:</p>
</div>
"""|>HTML

# ╔═╡ b5b430f4-6c0d-4ad0-b436-4b1ebea70ccc
begin
	err2=zeros(nr1-1,n)
	#h3= @. 1.0/(nRange1-1)
	plotBTCS_e=scatter(axis= :log,xlabel="\$ h\$",ylabel="successive Error in BTCS",legend= :bottomright, title="slution convergence")
    for i in 1:( nr1-1 )
        i1=index2.(res_BTCS[1].h_range,res_BTCS[i].h_range[2])
	    i2=index2.(res_BTCS[1].h_range,res_BTCS[i+1].h_range[2])
		for j in 2:nRange1[1]
	    err[i,j]= abs.(res_BTCS[i+1].all_time[i2[j],100]-res_BTCS[i].all_time[i1[j],100])
		end
    end
	lable2=["ξ=$ξ" for i in 1:1,ξ in res_BTCS[1].h_range[2:end] ]
	scatter(plotBTCS_e, h3[1:end-1], err[:,2:end], label=lable2 )
end

# ╔═╡ dda313d0-05b2-46f5-9948-3da07d886dcf
begin
	
	#nr1=length(nRange1)
    plotBTCS_t=scatter(xlabel="h",ylabel="FTCS times in end of fin ", title="slution for diffrent gride",legend= :bottomright)
    for i in 1:nr1-1 
	  scatter!(plotBTCS_t, res_BTCS[i].time_lapse, res_BTCS[i].all_time[end,:], label="n=$(nRange1[i])" )
    end
	plotBTCS_t
end

# ╔═╡ a030b0f5-6d4a-415c-99f5-dc65c56ef90b
"""<div style="direction: rtl"> <h2>   انجام محاسبات با روش: CTCS</h2><hr </div>
<div div style="direction: rtl;padding:14px">  
<p> روش CTCS به دلیل تفاوت در تجریه و معماری با سایر روش ها تفاوت جزئی دارد و این موضوع باعث شده تا نیاز به باز تعریف تعدادی از توابقع اولیه داشته باشیم (البته زیاد لزوم نداشت اما از این که توابع قبلی را جامع تر کنیم و تغییر سراسری بدیم راحت تر بود )   </p>
</div>

"""|>HTML

# ╔═╡ b2f8acf1-785f-49de-a4bd-0ce4642e8513
function initilize_variables_CTCS(n)
	return (
    Vector{Int64}(undef,3*(n-2)+4),
	Vector{Int64}(undef,3*(n-2)+4),
	Vector{Float64}(undef,3*(n-2)+4),
	zeros(Float64,n),
	zeros(Float64,n,n)
	)
end

# ╔═╡ 5b3fdf6b-11df-4662-aeca-5acd3aec7e33
function internal_coeffs_CTCS_n1!(I, J, V, C,P, n, h, τ, c1, k, c3, T_inf)
	for i= 2:(n-1)
		K=3*(i-2).+(1:3)# [3*n+1 ][3*o+1, 3*o+2, 3*o+3]
		I[K] .= i#[i, i, i]
		J[K] = [i-1, i, i+1]
		V[K] = [2*k*τ/c1/h^2, -2*(2*k*τ/c1/h^2+c3*τ/c1), 2*k*τ/c1/h^2]
		C[i] = 2*c3*τ*T_inf/c1
		P[i,i]=1
	end
end

# ╔═╡ 0c0888d7-c396-4770-b8fd-b666b8531eda
function right_bc_CTCS!(I, J, V, C, P, n, h, hc, k, T_inf)
        K= 3*(n - 2)+1 .+(1:3)
	    p=3*(n-3).+(1:3)
		I[K] .= n
		J[K] = [n-2, n - 1, n]
		V[K] =(k/(k+hc*h))*V[p]
	    C[n] = T_inf*hc*h/(k+hc*h)+C[n-1]
	    P[n,n-1]=(k/(k+hc*h))
	#return I, J, V, C
end

# ╔═╡ 44599fe2-877c-475d-ad0e-661cc2872b08
function left_bc_CTCS!(I, J, V, C, P, n)
        K= 3*(n-2) .+(1:1)
		#k=3*(i-2).+(1:3)
		I[K] .= 1
		J[K] = [1]
		V[K] = [1]
	    C[1] = 0
	    P[1,1]=0
	   # b[1] = T0
	#return I, J, V, C
end

# ╔═╡ 2e045543-272d-4508-a53c-af4d449680fc
bc_CTCS=(left_bc_CTCS!,right_bc_CTCS!)

# ╔═╡ f0812af2-77fa-45f2-8301-eb65eced109a
function generate_les_CTCS(n, hc, k, τ, c1, c3, T_inf,bc)
    h=dh(n)
	lbc,rbc=bc
	I, J, V, C, P= initilize_variables_CTCS(n)
	internal_coeffs_CTCS_n1!(I, J, V, C,P, n, h, τ, c1, k, c3, T_inf)
	lbc(I, J, V, C, P, n)
	rbc(I, J, V, C,P, n, h, hc, k, T_inf)
	return sparse(I, J, V) , C, P
end

# ╔═╡ f37dc87b-a04b-4a69-b2ca-fc5f7cc1a860
A_CTCS, C_CTCS, P_CTCS=generate_les_CTCS(n, hc, k, τ, c1, c3, T_inf,bc_CTCS)

# ╔═╡ 1422ce00-b96a-4967-9a54-c1330ff35043
function generate_all_time_CTCS(n, hc, k, τ, c1, c3,T0,T_inf,end_time,bc1,bc2)
	h_range=LinRange(0,1,n)
    time0=[T0;ones(n-1,1)*T_inf]
	A_CTCS, C_CTCS, P_CTCS = generate_les_CTCS(n, hc, k, τ, c1, c3, T_inf,bc1)
	
	time_lapse=collect(τ:τ:end_time)
	cloum=length(time_lapse)
	all_time=zeros(Float64,n,cloum)
	res_BTCS=generate_all_time_BTCS(n, hc, k, τ, c1, c3, T0, T_inf, 2*τ,bc2)
	all_time[:,2]=res_BTCS.all_time[:,1]
	all_time[:,1]=time0
	for i in (3:cloum)
		T2=A_CTCS*all_time[:,i-1]+P_CTCS*all_time[:,i-2]+C_CTCS
		all_time[:,i]=T2
		#time0=T
	end
     return (all_time=all_time , time_lapse=time_lapse, h_range=h_range)
end

# ╔═╡ d10602e0-16e9-42d6-a64e-b47a693a9705
all_time_CTCS=generate_all_time_CTCS(n, hc, k, τ, c1, c3,T0,T_inf,end_time,bc_CTCS,bc)

# ╔═╡ 7ecebffa-9a96-43a0-928d-d3b943879452
"""
<div style="direction: rtl"> <h3>   تحلیل خطا در روش CTCS</h3><hr </div>



<div div style="direction: rtl;padding:14px"> <p> همان طور که قابل مشاهده در نمودار است در این روش برای هیچ کدام از n های در نظر گرفته شده پایدار نبوده است.


   </p>



<p> </p>
</div>
"""|>HTML

# ╔═╡ 1f5246a6-e9db-4319-a9a2-ee0182007d43
 res_CTCS=[ generate_all_time_CTCS(ni, hc, k, τ, c1, c3,T0,T_inf,end_time,bc_CTCS,bc) for ni in nRange1 ]

# ╔═╡ ea73204c-42b8-43a7-a2b6-5e82fae660e1
begin
	
	#nr1=length(nRange1)
    plotB=scatter(xlabel="h",ylabel="CTCS in time=0.5", title="slution for diffrent gride")
    for i in 1:nr1-1 
	  scatter!(plotB, res_CTCS[i].h_range, res_CTCS[i].all_time[:,100], label="n=$(nRange1[i])" )
    end
	plotB
end

# ╔═╡ 880bc247-7440-41f0-89d1-28ea6874fae8
"""<div style="direction: rtl"> <h2>   انجام محاسبات با روش: کرنک نیکسون</h2><hr </div>


"""|>HTML

# ╔═╡ 4894234f-cff8-4bd4-8e6f-79ab43c81258
function initilize_variables_KN(n)
	return (
    Vector{Int64}(undef,3*(n-2)+3),
	Vector{Int64}(undef,3*(n-2)+3),
	Vector{Float64}(undef,3*(n-2)+3),
	Vector{Float64}(undef,3*(n-2)+3),
	zeros(Float64,n)
	)
end

# ╔═╡ 58c2cca9-88d8-4c9c-8205-b19d47ce29c3
function internal_coeffs_KN!(I, J, V_n, V_n2, C, n, h, τ, c1, k, c3, T_inf )
	for i= 2:(n-1)
		K=3*(i-2).+(1:3)# [3*n+1 ][3*o+1, 3*o+2, 3*o+3]
		I[K] .= i#[i, i, i]
		J[K] = [i-1, i, i+1]
		#V[K] = [2*k*τ/c1/h^2, -2*(2*k*τ/c1/h^2+c3*τ/c1), 2*k*τ/c1/h^2]
		V_n[K] = [k*τ, (-2*c1*h^2-2*k*τ-c3/2), k*τ]
		V_n2[K] = [-k*τ, (-2*c1*h^2+2*k*τ+c3/2), -k*τ]
		C[i] = -2*c3*τ*T_inf*h^2
	end
end

# ╔═╡ b81b5468-3191-4905-b50e-ba3832f452ae
function generate_les_KN(n, hc, k, τ, c1, c3, T_inf,bc)
    h=dh(n)
	lbc,rbc=bc
	I, J, V, C= initilize_variables(n)
	Cn=copy(C)
    Vn=copy(V)
	internal_coeffs_KN!(I, J, V, Vn, C, n, h, τ, c1, k, c3, T_inf )
	In =copy(I)
	Jn=copy(J)
	lbc(I, J, V, C, n)
	rbc(I, J, V, C, n, h, hc, k, T_inf)
	lbc(In, Jn, Vn, Cn, n)
	rbc(In, Jn, Vn, Cn, n, h, hc, k, T_inf)
	C[n]= C[n]-Cn[n]
	return sparse(I, J, V), sparse(In, Jn, Vn), C
end

# ╔═╡ 4938f673-c9c5-4aad-b358-10b0395e0036
A_kn1,A_kn2, C_kn=generate_les_KN(n, hc, k, τ, c1, c3, T_inf,bc)

# ╔═╡ d7a3dbac-736d-4ed8-b527-fa9ad314c7f4
T1=A_kn1\(A_kn2*(time0)+C_kn)

# ╔═╡ 1235041c-38d0-438e-a285-c15a70d10b3d
function generate_all_time_KN(n, hc, k, τ, c1, c3,T0,T_inf,end_time,bc)
    h_range=LinRange(0,1,n)
	time0=[T0;ones(n-1,1)*T_inf]
	A_kn1,A_kn2, C_kn = generate_les_KN(n, hc, k, τ, c1, c3,T_inf,bc)
	time_lapse=collect(τ:τ:end_time)
	cloum=length(time_lapse)
	all_time=zeros(Float64,n,cloum)
	for i in (1:cloum)
		T1=A_kn2\(A_kn1*(time0)+C_kn)
		all_time[:,i]=T1
		time0=T1
	end
	return (all_time=all_time , time_lapse=time_lapse, h_range=h_range)
end

# ╔═╡ 00a7750e-0618-43b8-9729-fc9c1393031d
res_KN= [ generate_all_time_KN(ni, hc, k, τ, c1, c3,T0,T_inf,end_time,bc) for ni in nRange1 ]

# ╔═╡ 90bcd55a-b47d-437f-8d5f-250c5b5b8f1a
slope(e,h)=?.log(e[2:end]/e[1:end-1])/log(h[2:end]/h[1:end-1])

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[compat]
Plots = "~1.26.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c9a6160317d1abe9c44b3beb367fd448117679ca"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.13.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "96b0bc6c52df76506efc8a441c6cf1adcb1babc4"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.42.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ae13fcbc7ab8f16b0856729b050ef0c446aa3492"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.4+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "9f836fb62492f4b0f0d3b06f55983f2704ed0883"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.0"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a6c850d77ad5118ad3be4bd188919ce97fffac47"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.0+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a6552bfeab40de157a297d84e03ade4b8177677f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.12"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "3f7cb7157ef860c637f3f4929c8ed5d9716933c6"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.7"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "737a5957f387b17e74d4ad2f440eb330b39a62c5"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "648107615c15d4e09f7eca16307bc821c1f718d8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.13+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "85b5da0fa43588c75bb1ff986493443f821c70b7"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "6f1b25e8ea06279b5689263cc538f51331d7ca17"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.1.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "23d109aad5d225e945c813c6ebef79104beda955"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.26.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "de893592a221142f3db370f48290e3a2ef39998f"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.4"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "c6c0f690d0cc7caddb74cef7aa847b824a16b256"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "995a812c6f7edea7527bb570f0ac39d0fb15663c"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.1"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "74fb527333e72ada2dd9ef77d98e4991fb185f04"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.1"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c3d8ba7f3fa0625b062b82853a7d5229cb728b6b"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.1"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "57617b34fa34f91d536eb265df67c2d4519b8b98"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.5"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─b0f3964c-9f15-4148-904f-815e0d78b017
# ╟─77b0c19c-0b5a-4b62-9b56-eb97b6434e2a
# ╟─7ce121eb-ff77-4c54-8d8f-07d286d603f0
# ╠═c529e753-f54a-4eb4-aa23-3ec975e27a84
# ╠═56eee02f-091a-4b5c-a850-08fb59f12aae
# ╠═9a1a0057-54a5-4c33-93ff-c3cdd70f1662
# ╠═eaa783cb-f687-4075-b320-3f7050636e04
# ╠═7e5e2c1a-6b51-43fd-83e4-78cd416615ff
# ╠═ac541cbb-c35c-490b-a0a6-4daef8cd3566
# ╟─aece47b6-9aa9-474d-ab82-da97eb53c2ad
# ╠═ee228b61-0538-49f0-ab30-74246458cf51
# ╟─5be587ad-7703-4430-a9eb-279a24d20fc0
# ╠═68b3dff3-1a99-4c31-ad3d-3b2455de69a6
# ╠═b90a98a0-e085-4f36-b415-c89a0bff1624
# ╠═5a4bc4e6-1bad-4a4a-bd73-00d56055d895
# ╟─60e0b122-25da-42a5-80e9-bf42bb6f6455
# ╠═b5eeb804-2cc8-436c-9acc-1277eb664195
# ╠═1b92dbcf-92c3-42e7-b7cb-c60bd5af3642
# ╠═59cd42b9-5792-429a-8d56-5e007d693150
# ╟─e4c74479-e565-4570-bd3a-fc4ec8c52b46
# ╟─b4f62d03-a975-4ce4-9de1-f6fd12944982
# ╠═0e70748b-8b49-44a2-9318-2938e727421d
# ╠═c2351b4f-89df-4a33-a073-92413245464c
# ╠═1e30de32-7004-4498-867c-7d8e028a4ee7
# ╟─2201cf8f-8e48-47db-a5f0-f6cc2265ef8f
# ╠═3b7e65df-07de-4c5c-8e42-e7159aeeb752
# ╟─1f3db9fb-0156-4f7d-bd6c-14fd4f512a45
# ╟─866d2f41-fdc3-4e49-ba91-49d95999eeed
# ╟─76164d36-37e4-4d6b-8d7b-b315c1c5c801
# ╠═be2591de-b5f1-4e16-b7d9-56c864c362db
# ╟─9bd7a2c0-5ab2-4781-bd25-6cf42edc0e96
# ╠═c5c42ec3-463a-4814-b9b3-eac861f1b7bc
# ╠═e9bfb747-0921-4cc7-ab6b-20ef741c3453
# ╠═b31236ce-3469-4759-9697-26ffa86566fc
# ╠═7f1cd37b-27fb-45fb-961e-9bfd2e979b55
# ╠═aea48b90-9784-46ec-83b8-e6757b428184
# ╠═fc63e35c-0b20-4a36-8760-543db4f0cac5
# ╟─c1cbc128-b47d-40b2-a969-ddfe919ca11c
# ╠═9298374c-c9af-4cef-ab61-e53922636d75
# ╠═a68b1cca-7642-4593-8ba0-04b36931c20a
# ╠═70ea54df-7d89-4411-90a4-56230db9c868
# ╟─3a950fd4-b631-4868-bcb7-07089d2eb4f9
# ╠═0b754307-067b-4ab8-b857-7b340c82eb50
# ╠═e2649731-93ce-4783-942f-82b6954bf061
# ╠═b2398175-5d89-4311-9651-0b494984a37a
# ╟─befb0951-9cc9-4f5c-ac6d-81ae08842b66
# ╠═fdc80205-b52b-4606-a6f9-16288c674299
# ╠═1878e5d0-ff15-4a6d-a2da-30f86bf9ab22
# ╟─e218cf86-2e2b-412f-b6ed-e5d6ba340f39
# ╠═f68d5d7e-3f64-4b04-9aab-31448a59fd3c
# ╟─f069dfeb-0cc6-49df-86c3-5aa60f807f99
# ╠═583a8f2c-913b-4eda-ad90-f66b3a38d87b
# ╟─7ff114c1-ebcd-48b8-8e5e-ec0e2168dfed
# ╠═cb06c89c-b0fc-43d4-b145-5724bfcbfb63
# ╟─3689c5df-404d-483d-91a1-c3df537deb62
# ╠═9e47327b-2e81-4b48-9331-09a5e2aba870
# ╠═c2a528ad-c689-4ce9-84dd-4fe532cb317d
# ╠═0a9f70c5-9be2-44d3-8541-7f0111977ca6
# ╠═c3da6280-8d6e-4e26-94d2-e1cd3a3e292b
# ╟─b7684934-c28f-4d12-bab5-7b521adfe9bb
# ╠═d1357278-e475-4e85-8491-306847558306
# ╠═8d192a41-7f8e-48ea-ba0c-7bd94872e1f8
# ╠═55b27c1d-8ae1-49b6-986a-eab21ea2383d
# ╠═71390fa2-e25b-4039-89ef-9f915e0b68d3
# ╠═8af39678-beaf-4d30-8d21-7704c8601219
# ╠═7880a698-2149-46bf-9010-b55bf32d22fb
# ╠═280e11f6-4edb-449d-a873-91ed7161dabb
# ╠═78808f18-5bc5-4ee4-825b-eadcf3c3ca45
# ╟─69fdbe90-7ce4-4a09-8b45-018d425d741b
# ╠═928a6042-3c76-47c3-8c76-7caf7397b981
# ╟─5aa48d02-15a0-4747-8dfe-c4cc76c86fbc
# ╠═b5b430f4-6c0d-4ad0-b436-4b1ebea70ccc
# ╠═dda313d0-05b2-46f5-9948-3da07d886dcf
# ╟─a030b0f5-6d4a-415c-99f5-dc65c56ef90b
# ╠═b2f8acf1-785f-49de-a4bd-0ce4642e8513
# ╠═5b3fdf6b-11df-4662-aeca-5acd3aec7e33
# ╠═0c0888d7-c396-4770-b8fd-b666b8531eda
# ╠═44599fe2-877c-475d-ad0e-661cc2872b08
# ╠═2e045543-272d-4508-a53c-af4d449680fc
# ╠═f0812af2-77fa-45f2-8301-eb65eced109a
# ╠═f37dc87b-a04b-4a69-b2ca-fc5f7cc1a860
# ╠═1422ce00-b96a-4967-9a54-c1330ff35043
# ╠═d10602e0-16e9-42d6-a64e-b47a693a9705
# ╟─7ecebffa-9a96-43a0-928d-d3b943879452
# ╠═1f5246a6-e9db-4319-a9a2-ee0182007d43
# ╠═ea73204c-42b8-43a7-a2b6-5e82fae660e1
# ╠═880bc247-7440-41f0-89d1-28ea6874fae8
# ╠═4894234f-cff8-4bd4-8e6f-79ab43c81258
# ╠═58c2cca9-88d8-4c9c-8205-b19d47ce29c3
# ╠═b81b5468-3191-4905-b50e-ba3832f452ae
# ╠═4938f673-c9c5-4aad-b358-10b0395e0036
# ╠═d7a3dbac-736d-4ed8-b527-fa9ad314c7f4
# ╠═1235041c-38d0-438e-a285-c15a70d10b3d
# ╠═00a7750e-0618-43b8-9729-fc9c1393031d
# ╠═90bcd55a-b47d-437f-8d5f-250c5b5b8f1a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
