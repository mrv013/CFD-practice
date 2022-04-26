### A Pluto.jl notebook ###
# v0.19.2

using Markdown
using InteractiveUtils

# ╔═╡ 9e25af28-7ba6-4db1-a659-36e5e43aeea1
using LinearAlgebra

# ╔═╡ 287d567a-7905-4801-a671-be58fcd5b8b0
using SparseArrays

# ╔═╡ 051ec1b2-33fe-4119-925a-ec1a34a20f63
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
 <p><b>
تکلیف سری سوم
</b>
مطالعه معادله موج و ضریب تقویت در 4 روش upwind, Lax, Lax-Wendroff, MacCormack
 <p>
 <p><b>عنوان درس:</b>
CFD
 <p>
 <p> <b>دانشجو:</b>
محمد رضا واعظی  <p>
 <p><b>استادر راهتما:</b>
دکتر نادران 
 <p>
 <p><b>تاریخ تحویل:</b>

1401/2/6
<p>
<a href="https://github.com/mrv013/CFD-practice.git">github repository</a>
<p>
(برای فعال شدن لینک صبر کنید)
<p>
<p>
https://github.com/mrv013/CFD-practice.git
<p>
<body></div>
"""|>HTML

# ╔═╡ 7ce121eb-ff77-4c54-8d8f-07d286d603f0
"""<div style="direction: rtl"> <h1>  مقدمه و پارامتر ها </h1></div>

<div  style="direction: rtl"> <p> در این تمرین قصد داریم تا مسئله و معادلات طرح شده را به چهار روش upwind, Lax, Lax-Wendroff, MacCormack تجزیه کرده و سپس معادلات را حل و در نهایت تحلیل خطا را انجام دهیم . </p>
<p>
معادله طرح شده به صورت زیر است :
</p>
</div>

"""|>HTML

# ╔═╡ 35e7ba37-23a2-4d03-8be0-b01997c71c55
md" 

```math
\begin{equation}
\frac{\partial u(x,t)}{\partial t}+
a \times \frac{\partial u(x,t)}{\partial x}=0
\end{equation}

```
"

# ╔═╡ e21545ec-f882-475e-ac99-2e4b3d6525c5
"""
<div  style="direction: rtl"> 
در دامنه از 0 تا L و زمان نامحدود با شرایط اولیه و مرزی در بخش مربوط به آن بحث شده است
</p>"""|>HTML

# ╔═╡ 0c71ef7b-dbaa-4918-a2ad-99cc49500757
"""
<div style="direction: rtl"> <p> برای انجام یک تحلیل CFD مناسب باید 3 بخش اصلی حفظ شود شامل 1- معادلات حاکم ، 2- تجزیه معادلات ، 3- تحلیل خطا 
</p>
<p>
سایر بخش ها مانند شبکه محاسباتی  با توجه به سادگی مسئله به صورت مختصر مورد توجه قرار گرفته اند در ادامه در هر فصل و زیر فصل توضیح داده شده است 
</p>
</div>

"""|>HTML

# ╔═╡ fe026e1a-8bee-46cb-8f0a-f9e390321f23
"""
<div style="direction: rtl"> <p>حال در ابتدا پارامتر های مسئله را تعریف میکنیم همان طور که  از فرمول مشخص است ما تنها یک متغیر a در فرمول داریم. </p> </div>

"""|>HTML

# ╔═╡ 56eee02f-091a-4b5c-a850-08fb59f12aae
 GC.gc()

# ╔═╡ 9a1a0057-54a5-4c33-93ff-c3cdd70f1662
a=0.5

# ╔═╡ aece47b6-9aa9-474d-ab82-da97eb53c2ad
"""

<div style="direction: rtl;padding:14px";
 #text-align:center;> <p>τ , h به ترتیب فاصله مکانی و زمانی مسئله هتند و در این بخش به تعریف آن ها اقدام میکنیم  </p> </div>

"""|>HTML

# ╔═╡ ee228b61-0538-49f0-ab30-74246458cf51
τ=0.01

# ╔═╡ 5be587ad-7703-4430-a9eb-279a24d20fc0
"""
<div div style="direction: rtl;padding:14px"> <p> برای پیدا کردن h را وابسته به تعداد نقاط مورد نظر برای حل مسئله را می کنیم سپس با تابع تعریف شده مقدار h را تعیین میکنیم  </p> 
<p> (
Δx همان h است که در برخی از توابع در از آن استفاده شئده است
)</p>

</div>

"""|>HTML

# ╔═╡ 68b3dff3-1a99-4c31-ad3d-3b2455de69a6
n=11

# ╔═╡ c34595ea-f229-47e8-aa01-b187bc9f90c2
L=1

# ╔═╡ b90a98a0-e085-4f36-b415-c89a0bff1624
dh(n,L)=L/(n-1)

# ╔═╡ 5a4bc4e6-1bad-4a4a-bd73-00d56055d895
h=dh(n,L)

# ╔═╡ 94915014-48b1-4c3c-93bb-5f467c3d3502
Δx=L/(n-1)

# ╔═╡ 28fb4806-a658-4489-9296-c36b1e1c504e
md"
```math
\begin{equation}
ν= \frac{τ.a}{h}
\end{equation}
```
"

# ╔═╡ 118ad3b0-fec2-4b86-9e58-c1f099602819
ν=a*τ/h

# ╔═╡ 60e0b122-25da-42a5-80e9-bf42bb6f6455
"""
<div div style="direction: rtl;padding:14px"> <p>در پایان زمان پایانی برای حل مسئله را مشخص و مقادیر اولیه دما و دمای محیط را مشخص میکنیم </p> </div>

"""|>HTML

# ╔═╡ 39d600fd-a18d-4340-96af-2b2ccd3f50fc
end_time=5

# ╔═╡ e4c74479-e565-4570-bd3a-fc4ec8c52b46
"""<div style="direction: rtl"> <h1>  حل عددی </h1></div>


"""|>HTML

# ╔═╡ c2df185f-5284-41e9-a662-8246b27d684b
"""<span style="direction: rtl"> <h2> شرایط مرزی و شرایط اولیه </h2></span>

"""|>HTML

# ╔═╡ 509b28b6-5152-4dc3-b2b8-4ec4a2101df8
"""
<div div style="direction: rtl;padding:14px"> <p> شرایط مرزی در تمام روش های محاسباتی یکسان هستند به همین علت توابع آن ها قبل از هر چیز تعریف شده است  </p> </div>

"""|>HTML

# ╔═╡ 9819c30a-08c8-4ca4-9f4b-e158b35a89d8
"""
<div div style="direction: rtl;padding:14px">  
<p> 
<b>
شرایط اولیه
</b>
باید طوری تعریف شود تا در هر نقطه از شبکه مقدار مشخصی نسبت داده شود یعنی بردار شرایط اویه باید به شکل شبکه وابسته باشد همچنین شبکه میتواند طوری باشد که شرایط اولیه به درستی دیده نشود مثلا اگر n=5 باشد شرایط اولیه حالت 3 بدرستی منعکس نمسشود و عملا مقدار 1 در آن دیده نمشود این محدودیت نه از سر استنسیل که از سر شکل شبکه و شرایط اولیه است.</p>

</p>  </p>

<p>  تعریف ماتریس شرایط اولیه(در زمان صفر):</p>
</div>
"""|>HTML

# ╔═╡ cfbc4bb5-56bb-4659-9036-71e7e092744f
md" 
```math
\begin{equation}
1:u(x,0)=
\left\{ 
\begin{aligned}
  1 \qquad 0\leq x<0.25\\
  0 \qquad 0.25\leq x\leq 1
\end{aligned}
\right.
\end{equation}
```
```math
\begin{equation}
2:u(x,0)=sin(4πx)
\end{equation}
```
```math
\begin{equation}
3:u(x,0)=
\left\{ 
\begin{aligned}
  0 \qquad 0\leq x<0.2\\
  1 \qquad 0.2\leq x<0.3\\
  0 \qquad 0.3\leq x\leq 1
\end{aligned}
\right.
\end{equation}

```


"

# ╔═╡ 14146054-6afb-456c-8db9-94fcde00d331
index(ξ,Δξ)=floor(Int, ξ/Δξ)+1

# ╔═╡ 1c68b2a3-62c6-4449-9e2c-6c69b9a3df04
"""
<div div style="direction: rtl;padding:14px">  
</p>  </p>

<p>  تعریف ماتریس شرایط اولیه(در زمان صفر) برای هر 3 نوع شرایط اولیه خواسته شده برای به عهده تابعی تعریف شده قرار میدهیم به این صورت نیاز به تغییر و تعریف شرایط اولیه به صورت جدا نداریم . ty  که یکی از ورودی های تابع inatioal_value  هست نوع شرایط اولیه خواسته شده را تعریف میکند </p>
</div>
"""|>HTML

# ╔═╡ ec51d580-7305-4be4-9e64-6507c1230cac
function inational_value!(L,n,ty)
	time0=zeros(n,1)
	Δx=L/(n-1)
	if (ty==1)
	  for i=1:n
		  x=(i-1)*Δx
		 if x<0.25 time0[i]=1 end
	  end
	elseif (ty==3)
	  for i=1:n
		 x=(i-1)*Δx
		 if x>=0.2 && x<0.3 time0[i]=1 end
	  end
	elseif (ty==2)
	  for i=1:n
		 x=(i-1)*Δx
		 time0[i]=sin(4*pi*x)
	  end
	end
return time0
	end		

# ╔═╡ eb548fea-9c56-4b39-9d2a-2caa9926b80f
u0_1=inational_value!(L,n,1)

# ╔═╡ dd23acef-14fe-41c3-b991-b61256b24d75
u0_2=inational_value!(L,n,2)

# ╔═╡ 1c3a2a15-a995-4124-9bf2-c61b5c62af23
u0_3=inational_value!(L,n,3)

# ╔═╡ 4453c8b4-e94c-4095-a3cf-4ac64be7d35d
"""
<div div style="direction: rtl;padding:14px"> 

<p>
<b>
 شرایط مرزی:
</b>
</p>

<p> با توجه به نوع مسئله (a>0) نیازی به شرایط مرزی در سمت راست)(x=L) معادله نداریم اما در سمت راست به دلیل جلو رفتن موج نیاز به مقادیر جدید داریم به همین دلیل مسئله عملا با یک شرط مرزی سمت چپ قابل حل خواهد بود. اما باید به این نکته توجه کرد که در استنسیل هایی که یک گام مکانی جلو تر در ان ها دیده میشود در x=L عبارتی برای نقطه L+Δx در ان ها هست که باید به نحوی از معادلات حذف  و جاگذاری شود . </p>
<p>  ما در این حل شرایط مرزی در 
<b>
سمت چپ
</b>
را برابر و ثابت  (u(x,0)= u<sub>0</sub> )   فرض کرده ایم. برای شرایط اولیه حالت 3 میتوان شرایط مرزی را متغیر فرض شود اما در آن صورت انتقال موج مناست در آن دیده نمشود. مقدار ثابت فرض شده میتواند هر مقداری داشته باشد و ما در این حل آن را برابر مقدار در شرایط اولیه در نقطه x=0 فرض کرده ایم این یعنی برای هر شرایط اولیه مقدار ثابت تغییر میکند اما ماتریس ظرایب ما ثابت میماند.</p>
<p>
 تابع تعریف ضرایب در نقاط مرزی:
</p>

</div>


"""|>HTML

# ╔═╡ bfc059e1-2561-46d9-b2c1-71adbdbdbc7e
md" 


```math
\begin{equation}
u_{0}^{n}= u_{0}^{n+1}=0
\end{equation}
```
"

# ╔═╡ 022e31c3-8864-4106-9aa4-6ebd49124814
function left_bc!(I, J, V, n)
        K= 3*(n-2) .+(1:1)
		#k=3*(i-2).+(1:3)
		I[K] .= 1
		J[K] = [1]
		V[K] = [1]
	   # C[1] = 0
	   # b[1] = T0
	#return I, J, V, C
end

# ╔═╡ b5556da7-056d-493f-bded-40d73376bd58
"""
<div div style="direction: rtl;padding:14px"> <p> همان طور که گفته شد  تابع در
<b>
سمت راست
</b>
به شرایطی برای حل نیاز ندارد مثلا در روش upwind کلا نیازی به شرایط مرزی سمت راست ندارد اما در سایر موارد(بدلیل وجود گام زمانی بعدی در پیکره محاسباتی) نیاز به نوعی شرط مرزی داریم. ما برای این حالت برای نقطه L با استفاده از روش upwind اقدام به تجزیه کرده کرده ایم تا نقطه L+Δx در شرایط مرزی سمت راست در پیکره محاسباتی دخیل نشود. </p>
<div div style="direction: rtl;padding:14px"> 
<p>  شرایط مرزی در
<b>
سمت راست
</b>
 (نقطه انتهایی شبکه) تعریف شده بر اساس روش upwind و معادله اصلی در نقطه L: </p>
</div>


</div>
"""|>HTML

# ╔═╡ b3efc4cc-25a0-46a7-b1e4-7d6d15194f32
md" 

```math
\begin{equation}
\frac{\partial u(x,t)}{\partial t}+
a \times \frac{\partial u(x,t)}{\partial x}=0
\end{equation}

```
"

# ╔═╡ ae7ba8c0-d12c-45c2-99f4-863b52f314b7
"""
<div div style="direction: rtl;padding:14px"> <p>  تجزیه معادله  شرط مرزی سمت راست و تعریف تابع ضرایب در نقات مرزی سمت راست :</p>
</div>


</div>
"""|>HTML

# ╔═╡ 1017451b-46c3-4009-8fd8-ba257db39ac2
md" 

```math
\begin{equation}
u_{L}^{ n+1}= \nu\times u_{L-1}^{n}+
(1-\nu)\times u_{L}^{ n}
\end{equation}

```
"

# ╔═╡ a2daaecd-b154-4363-9feb-2c87d3300538
function right_bc!(I, J, V, n, ν)
        K= 3*(n - 2)+1 .+(1:2)
	    #p=3*(n-3)
		I[K] .= n
		J[K] = [n - 1, n]
	    V[K] = [ν, 1-ν]
	    #V[K] = V[p+1:p+2]+ V[p+3]*[ν, 1-ν]
	    #V[K] = [k*τ/c1/h^2, (1.0 -2*k*τ/c1/(h^2)-c3*τ/c1)]+k*τ/c1/(h^2)*[-2*hc*h/k 1]
		#V[K] =(k/(k+hc*h))*V[p]
	    #C[n] = (V[p+3])*2*T_inf*hc*h/k + C[n-1]
	return I, J, V
end

# ╔═╡ 65de7f54-9384-4124-b353-fa1ef495a1a6
bc = (left_bc!,right_bc!)

# ╔═╡ 9c060b1d-a057-489d-844f-69992df65bc1


# ╔═╡ fef9cd8b-05a2-403e-a092-81d9960feeff
"""<span style="direction: rtl"> <h2> تعریف کتاب خانه ها و حافظه  </h2></span>

"""|>HTML

# ╔═╡ 503388cf-eeb1-4c42-b04f-0263bcfb316e
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
	#zeros(Float64,n),
	#zeros(Float64,n)
	)
end

# ╔═╡ 1f3db9fb-0156-4f7d-bd6c-14fd4f512a45
"""
<div div style="direction: rtl;padding:14px"> <p> در ادامه برای هر چهار روش MacCormack, Lax-Wendroff, Lax, upwind جدا گانه توابع و بررسی خطا وابسته به بازه های زمانی و مکانی را انجام داده و بررسی آن ها را تکمیل میکنیم. تمام این روش ها صریح هستند و پیچیدگی و حجم محاسباتی روش های ضمنی را نداریند اما باید رد باره بازه ها در آن ها دقت به خرج داد .  </p> </div>

"""|>HTML

# ╔═╡ 866d2f41-fdc3-4e49-ba91-49d95999eeed
"""<span style="direction: rtl"> <h2> انجام محاسبات با روش upwind</h2></span>


"""|>HTML

# ╔═╡ 9d190dbf-afcf-4293-97c5-13d2d5b7e5e8
"""
<div div style="direction: rtl;padding:14px"> <p>    </p>
</div>

"""|>HTML

# ╔═╡ c76ae0c1-ed35-4547-9fef-8f3c8a464092
"""<span style="direction: rtl"> <h3> تجزیه معادلات حاکم به روش upwind</h3></span>

"""|>HTML

# ╔═╡ 6743e0ca-b724-4822-94df-e8254853c842
"""
<div div style="direction: rtl;padding:14px"> <p>  معادله زیر تجزیه به روش upwind است: </p>
</div>

"""|>HTML

# ╔═╡ 9b33684a-337d-424b-a632-881fbe4e7571
md" 

```math
\begin{equation}
u_{i}^{ n+1}= \nu\times u_{i-1}^{n}+
(1-\nu)\times u_{i}^{ n}
\end{equation}

```
"

# ╔═╡ 007f9a43-f4c8-4d3e-96f9-895de84efa97
md" 

```math
\begin{equation}
G_{m}= 1+ \nu\times (cos(\beta)-1)+
i\nu\times sin(\beta)

\end{equation}

```
"

# ╔═╡ a56a5e2f-4302-438a-8921-c31e13674b22
begin
G_m_upwind(a,b)=1+a*(cos(b)-1)+(a*sin(b))im
b=range(0, stop=π, length=100)
plo_G_m_upwind=scatter( title="polar plot G_m upwind",legend= :bottomright, proj=:polar, rlims=(0,1.5),reuse=true)
for a=1:-0.25:0.25
    #println(a)
    #θ=nothing
    #G_m=nothing
    θ=@.angle(G_m_upwind(a,b))
    G_m=@.abs(G_m_upwind(a,b))
    plot!(plo_G_m_upwind,θ, G_m,label="ν=$a")
end
	
plot(plo_G_m_upwind)
end

# ╔═╡ 6d92754d-c4f0-4b34-ad0b-9f1b57391331


# ╔═╡ 8fc18020-e7bf-4060-abcb-5a4ff3c636df
"""<span style="direction: rtl"> <h3> کد ها و مراحل محاسبات کامپیوتری روش upwind  </h3></span>

"""|>HTML

# ╔═╡ 76164d36-37e4-4d6b-8d7b-b315c1c5c801
"""
<div div style="direction: rtl;padding:14px"> <p>  هر روش محاسبه 4 نوع تابع نیاز دارد اول باید ماتریس ظرایب را به روش OCC تعریف میکنیم سپس تابعی برای معرفی شرایط مرزی تعریف میکنیم و بعد به سراغ دو تابع بکب برای انجام محاسبات در دامنه مکان با یک قدم در زمان و دیگری برای گسترش محاسبات تا end _time  استفاده میکنیم    </p>
<p> تابع تعریف ضرایب برای نقات میانی :</p>
</div>

"""|>HTML

# ╔═╡ be2591de-b5f1-4e16-b7d9-56c864c362db
function internal_coeffs_upwind!(I, J, V, n, ν)
	for i= 2:(n-1)
		K=3*(i-2).+(1:3)# [3*n+1 ][3*o+1, 3*o+2, 3*o+3]
		I[K] .= i#[i, i, i]
		J[K] = [i-1, i,i+1]
		V[K] = [ν, 1-ν,0]
		#C[i] = c3*τ*T_inf/c1
	end
	#return I, J, V#, C
end

# ╔═╡ 7f1cd37b-27fb-45fb-961e-9bfd2e979b55
"""
<div div style="direction: rtl;padding:14px">  
<p> تعریف تابع برای محاسبه تابع Dense:</p>
</div>
"""|>HTML

# ╔═╡ fc63e35c-0b20-4a36-8760-543db4f0cac5
function generate_les_upwind(n, ν, bc)
    h=dh(L,n)
	lbc,rbc=bc
	I, J, V= initilize_variables(n)
	internal_coeffs_upwind!(I, J, V, n, ν)
	lbc(I, J, V, n)
	rbc(I, J, V, n, ν)
	return sparse(I, J, V) #, C
end

# ╔═╡ a68b1cca-7642-4593-8ba0-04b36931c20a
A_upwind = generate_les_upwind(n, ν, bc)

# ╔═╡ 70ea54df-7d89-4411-90a4-56230db9c868
det(A_upwind)

# ╔═╡ 3a950fd4-b631-4868-bcb7-07089d2eb4f9
"""
<div div style="direction: rtl;padding:14px">  
<p>  تعریف تابع برای محاسبه در کل بازه زمانی (همان طور که مشخص است یکی از ورودی های این تابع end_time  است ):</p>
</div>
"""|>HTML

# ╔═╡ 0b754307-067b-4ab8-b857-7b340c82eb50
function generate_all_time_upwind(n, u0, L, a, ν, end_time, bc)
	h_range=LinRange(0,L,n)
	τ=ν*h_range[2]/a
    time0=u0
	A = generate_les_upwind(n, ν, bc)
	time_lapse=collect(τ:τ:end_time)
	cloum=length(time_lapse)
	all_time=zeros(Float64,n,cloum)
	for i in (1:cloum)
		u1=A*(time0)
		all_time[:,i]=u1
		time0=u1
	end
	return (all_time=all_time , time_lapse=time_lapse, h_range=h_range)
end

# ╔═╡ e2649731-93ce-4783-942f-82b6954bf061
u1_upwind=generate_all_time_upwind(n,inational_value!(L,n,1),L, a, ν,end_time,bc)

# ╔═╡ b2398175-5d89-4311-9651-0b494984a37a
plo1s=scatter(u1_upwind.h_range,u1_upwind.all_time[:,2],xlabel="ξ",ylabel="θ", label="n=$n")

# ╔═╡ befb0951-9cc9-4f5c-ac6d-81ae08842b66
"""
<div style="direction: rtl"> <h3>   تحلیل خطا در روش upwind</h3><hr </div>
<p>
<b>
ν ثابت، n متغیر
</b>
</p>

<div div style="direction: rtl;padding:14px"> <p> در این بخش برای تحیلی خطا به 3 نوع عمل میکنیم
<b>نوع اول: </b>
با استفاده از توابع تعریف شده مسئله را برای چنداز n (تعداد نقاط در طول فین)با مقدار ثابت از tau انجام  و مقایسه مکنیم. 
<b>نوع دوم: </b>
خطای sucsecive  را برای نوع اول در نقاط مشخص محاسبه میکنیم.
<b>نوع سوم: </b>
  محاسبات را برای چند بازه زمانی(τ) و با تغییر   n  انجام میدهیم و باهم مقایسه میکنیم. در پایان نوع 2 و 3 را برای
<b>
n ثابت، ν متغیر 
</b>
اجرا میکنیم     
   </p>



<p> </p>
</div>
"""|>HTML

# ╔═╡ fdc80205-b52b-4606-a6f9-16288c674299
nRange1=[1,2,4,8].*(n-1).+1

# ╔═╡ 1878e5d0-ff15-4a6d-a2da-30f86bf9ab22
length(nRange1)

# ╔═╡ 2b7a7b49-97f8-4c02-aa0f-b855bc90b268
"""
<div div style="direction: rtl;padding:14px">  
<p><b> شرایط اولیه نوع دوم  (ν ثابت، n متغیر) </b>  </p>
<p>خطا برای شرایط اولیه نوع  دوم:  </p>
</div>
"""|>HTML

# ╔═╡ e218cf86-2e2b-412f-b6ed-e5d6ba340f39
"""
<div div style="direction: rtl;padding:14px">  
<p> محاسبات برای n موجود در مجموعه nRange:</p>
</div>
"""|>HTML

# ╔═╡ f68d5d7e-3f64-4b04-9aab-31448a59fd3c
res_upwind= [generate_all_time_upwind(ni,inational_value!(L,ni,2),L, a, ν,end_time,bc) for ni in nRange1 ]

# ╔═╡ f069dfeb-0cc6-49df-86c3-5aa60f807f99
"""
<div div style="direction: rtl;padding:14px">  
<p> نمایش خطای نوع اول برای n موجود در مجموعه nRange:</p>
</div>
"""|>HTML

# ╔═╡ 583a8f2c-913b-4eda-ad90-f66b3a38d87b
begin
	time_step_want=20
	nr1=length(nRange1)
    plotc=scatter(xlabel="h",ylabel="upwind in time=0.5", title="slution for diffrent gride")
    for i in 1:nr1 
	  scatter!(plotc, res_upwind[i].h_range, res_upwind[i].all_time[:,time_step_want], label="n=$(nRange1[i])" )
    end
	plotc
end

# ╔═╡ 3689c5df-404d-483d-91a1-c3df537deb62
"""
<div div style="direction: rtl;padding:14px">  
<p> در این مرحله برای محاسبه خطای successive ابتدا باید نقاط ثابتی(ما نقاط مکانی در n=11 را انتخاب کردیم) را در طول فین انتخاب کنیم و شماره نقاط را در هر n (n هارا متفاوت در مجموعه nReng1) پیدا کرده و سپس به محاسبه خطا و رسم نمودار پرداخت . تابع تعریف شده index  به این منظور است .  </p>
</div>
"""|>HTML

# ╔═╡ 2d92c93f-4e09-48a3-b92a-fb6bc9d49d24
slope(e,h)=@. log(e[2:end]/e[1:end-1])/log(h[2:end]/h[1:end-1])

# ╔═╡ 561a0b16-6865-4381-bb48-9dfa79e35c37
"""
<div div style="direction: rtl;padding:14px">  
<p> کاهش h باعث میشود که هارمونیک های بیشتری دیه شود در واقع در نمودار Gm فاز های بیشتری دیده میشوند و استهلاک آن ها نیز اثر میکند و باعث عدم همگرایی میشود. از طرفیبا کاهش h و ثابت ماندن ν گام زمانی افزایش میابد. با افزایش بیش از حد گام زمانی به روز رسانی نقاط از سرعت موج جلو میزند  با عث خطا بیشتر میشود .  </p>
</div>
"""|>HTML

# ╔═╡ 6f2f21c1-e9a7-4db4-a1c2-6204abe2ae48
"""
<div div style="direction: rtl;padding:14px">  
<p> در این مرحله با توجه به کد نمودار ها تابع ایی تعریف میکنیم تا بتوان روش های دیگر را راحت تر به تصویر کشید.   </p>
</div>
"""|>HTML

# ╔═╡ e5659b9e-5d9a-43ed-8a67-3b2245ffb138
time_want=1

# ╔═╡ 0e5f9639-c5a6-4a17-b077-8ad37361de80
"""
<div div style="direction: rtl;padding:14px">  
<p><b> شرایط اولیه نوع اول  (ν ثابت، n متغیر) </b>  </p>
<p>خطا برای شرایط اولیه نوع  اول :  </p>
</div>
"""|>HTML

# ╔═╡ ba0a7ab7-c687-4a87-9619-2e3e1730199a
res_upwind_err_1= [generate_all_time_upwind(ni,inational_value!(L,ni,1),L, a, ν,end_time,bc) for ni in nRange1 ]

# ╔═╡ ef11bf7a-33b6-4ef1-afa6-3077e5efd9ac


# ╔═╡ 47ac2a22-20c7-4a43-95ae-f9f7195a71a1
"""
<div div style="direction: rtl;padding:14px">  
<p><b> شرایط اولیه نوع سوم (ν ثابت، n متغیر)  </b>  </p>
<p>خطا برای شرایط اولیه نوع  سوم:  </p>
</div>
"""|>HTML

# ╔═╡ 07cb20d4-1571-4bec-9dae-00017e03ce07
res_upwind_err_3= [generate_all_time_upwind(ni,inational_value!(L,ni,3),L, a, ν,end_time,bc) for ni in nRange1 ]

# ╔═╡ 94dce462-9662-402f-85d3-074df56d96ab
"""
<div div style="direction: rtl;padding:14px">  
<p> 
<b>
n ثابت، ν متغیر 
</b>
</p>
<p> 
در این بخش با ثابت نگه داشتن ν و تغییر n در واقع τ را تغییر میدهیم.این نمودار های با معنی تر از نمودار های قبلی هستند و بهتر تفصیر میشوند.
</p>
</div>
"""|>HTML



# ╔═╡ 8c85bfb2-d5bd-45d5-8368-ab25f94ebc45
ν_Range=[1,0.5, 0.25, 0.125]*ν

# ╔═╡ 4699894f-4560-4658-9b4a-960eb7fb6b5f
method_name= "upwind"

# ╔═╡ 8a790d36-de2f-485c-b476-fee78772c796
time=[0.1, 0.2, 0.3, 0.4, 0.5]

# ╔═╡ 3d874be8-61c8-4620-bd2b-fbdb00a9f3ab
index_tim(ξ,Δξ)=floor(Int, ξ/Δξ)

# ╔═╡ 0a9f70c5-9be2-44d3-8541-7f0111977ca6
begin
	err=zeros(nr1-1,n)
	h3= @. 1.0/(nRange1-1) #[0.1, 0.05, 0.025, 0.0125]
	plo3=scatter(axis= :log,xlabel="\$ h\$",ylabel="successive Error upwind", title="slution convergence, time=$(time_step_want*τ) s",legend= :bottomright)
    for i in 1:( nr1-1 )
		time_step_want=index_tim.(time_want,res_upwind[i].time_lapse[1])
        i1=index.(res_upwind[1].h_range,res_upwind[i].h_range[2])
	    i2=index.(res_upwind[1].h_range,res_upwind[i+1].h_range[2])
		for j in 2:nRange1[1]
	    err[i,j]= abs.(res_upwind[i+1].all_time[i2[j],time_step_want]- res_upwind[i].all_time[i1[j],time_step_want])
		end
    end
	lable=["x=$ξ" for i in 1:1,ξ in res_upwind[1].h_range[2:end]]
	scatter(plo3, h3[1:end-1], err[:,2:end], label=lable )
end


# ╔═╡ ec900552-34c8-4486-b37c-a32536f61cfb
begin
	ploupwind_slope=plot(xlabel="\$ h\$",ylabel="Slope",title="slution convergence upwind, time=$(time_step_want*τ) s",legend= :topleft)
    for i in 2:nRange1[1]#length(nRange1)
		plot!(ploupwind_slope,h3[1:end-2], slope(err[:,i],h3[1:end-1]),label="x=$(res_upwind[1].h_range[i])",lw=2)
    end
	ploupwind_slope
end

# ╔═╡ ef5ddb37-91fd-446b-8333-e50a48008fc2
function plot_x_seris(ν, res_upwind, time_want, nRange1,method_name)
	#time_step_want=20
	nr1=length(nRange1)
    plotc=scatter(xlabel="h",ylabel="u", title="slution for diffrent 
    gride $method_name , time=$(time_want) s, ν=$ν")
    for i in 1:nr1
	  time_step_want=index_tim.(time_want,res_upwind[i].time_lapse[1])
	  scatter!(plotc, res_upwind[i].h_range, 
      res_upwind[i].all_time[:,time_step_want], label="n=$(nRange1[i])" )
    end
	#plotc
	
    err=zeros(nr1-1,n)
	h3= @. 1.0/(nRange1-1) #[0.1, 0.05, 0.025, 0.0125] 
	#τ=res_upwind[1].time_lapse[2]
	plo3=scatter(axis= :log,xlabel="\$ h\$",ylabel="successive Error", title="solution convergence $method_name, time=$(time_want) s, ν=$ν",legend= :bottomright)
    for i in 1:( nr1-1 )
		time_step_want=index_tim.(time_want,res_upwind[i].time_lapse[1])
        i1=index.(res_upwind[1].h_range,res_upwind[i].h_range[2])
	    i2=index.(res_upwind[1].h_range,res_upwind[i+1].h_range[2])
		for j in 2:nRange1[1]
	    err[i,j]= abs.(res_upwind[i+1].all_time[i2[j],time_step_want]-res_upwind[i].all_time[i1[j],time_step_want])
		end
    end
	lable=["x=$ξ" for i in 1:1,ξ in res_upwind[1].h_range[2:end] ]
	successive_Error=scatter(plo3, h3[1:end-1], err[:,2:end], label=lable )

	plot_slope=plot(xlabel="\$ h\$",ylabel="Slope",title="solution convergence $method_name, time=$(time_want) s, ν=$ν",legend= :topleft)
    for i in 2:nRange1[1]#length(nRange1)
		
		plot!(plot_slope,h3[1:end-2], slope(err[:,i],h3[1:end-1]),label="x=$(res_upwind[1].h_range[i])",lw=2)
    end
	
	
    return (solution=plotc, successive_Error=successive_Error, plot_slope=plot_slope)
end

# ╔═╡ 52e072cd-a21e-4a83-88b8-3fe3678c0fd8
plot_x_seris(ν, res_upwind, time_want, nRange1,method_name)

# ╔═╡ 467bf194-1daa-44d4-a1c7-78863d293555
plot_x_seris(ν, res_upwind_err_1, time_want, nRange1,method_name)

# ╔═╡ 5e759c6c-2ce9-4c1d-89de-01a700b015de
plot_x_seris(ν, res_upwind_err_3, time_want, nRange1,method_name)

# ╔═╡ 2bb874d6-2018-4c07-930f-d040bcc5eba7
"""
<div div style="direction: rtl;padding:14px">  
<p><b> شرایط اولیه نوع دوم (n ثابت، ν متغیر ) </b>  </p>
<p>خطا برای شرایط اولیه نوع  دوم:  </p>
</div>
"""|>HTML

# ╔═╡ ae7a2f50-ae54-49fa-8119-0ae062ee1a30
res_upwind_t= [ generate_all_time_upwind(n,inational_value!(L,n,2),L, a, νi,end_time,bc) for νi in ν_Range ]

# ╔═╡ dc43d2af-2975-4c94-80aa-c0cf7c189554
x_step_t=5

# ╔═╡ 201ff5bc-8d24-4400-bcbc-d16262c22347
begin
	n_time=length(time)
	n_ν=length(ν_Range)
	err_t=zeros(n_ν-1,n_time)
	#h3= @. 1.0/(nRange1-1)
	plot_upwind_t=scatter(axis= :log,xlabel="ν",ylabel="successive Error",legend= :bottomright,title="solution convergence upwind, x=$((x_step_t-1)*1/(nRange1[1]-1))*L ")
    for i in 1:( n_ν-1 )
        i1=index_tim.(time,ν_Range[i])
	    i2=index_tim.(time,ν_Range[i+1])
		for j in 1:n_time
	    err_t[i,j]= abs.(res_upwind_t[i+1].all_time[x_step_t,i2[j]]-res_upwind_t[i].all_time[x_step_t,i1[j]])
		end
    end
	lable_t=["t=$ξ" for i in 1:1,ξ in time ]
	scatter(plot_upwind_t, ν_Range[1:end-1], err_t[:,2:end], label=lable_t )
end

# ╔═╡ 061d2b74-96e1-43b8-b252-13911ba9639a
"""
<div div style="direction: rtl;padding:14px">  
<p> همان طور که در نمودار شکل زیر قابل مشاهده است با کاهش ν خطای تمام نقاط به هم همگرا میشود.  </p>
</div>
"""|>HTML

# ╔═╡ 89aea175-f4e1-4777-ba28-0c684a959295
begin
	ploupwind_slope_t=plot(xlabel="ν",ylabel="Slope",title="slution convergence upwind, x=$((x_step_t-1)*1/(nRange1[1]-1))*L ",legend= :bottomright)
    for i in 2:n_time
		plot!(ploupwind_slope_t,ν_Range[1:end-2], slope(err_t[:,i],ν_Range[1:end-1]),label="t=$(time[i])",lw=2)
    end
	ploupwind_slope_t
end

# ╔═╡ 27883575-5235-40b1-97aa-d1779b5f35b6
"""
<div div style="direction: rtl;padding:14px">  
<p> در این مرحله با توجه به کد نمودار ها تابع ایی تعریف میکنیم تا بتوان روش های دیگر را راحت تر به تصویر کشید.   </p>
</div>
"""|>HTML

# ╔═╡ 76fac5ab-9676-4bfe-a08b-4ce75ea98774
function plot_ν_seris(time, ν_Range, res_upwind_t, x_step_t, nRange1,method_name)
	n_time=length(time)
	n_ν=length(ν_Range)
	err_t=zeros(n_ν-1,n_time)
	#h3= @. 1.0/(nRange1-1)
	plot_upwind_t=scatter(axis= :log,xlabel="ν",ylabel="successive Error",title="slution convergence $method_name, x=$((x_step_t-1)*1/(nRange1[1]-1))*L ",legend= :bottomright)
    for i in 1:( n_ν-1 )
        i1=index_tim.(time,ν_Range[i])
	    i2=index_tim.(time,ν_Range[i+1])
		for j in 1:n_time
	    err_t[i,j]= abs.(res_upwind_t[i+1].all_time[x_step_t,i2[j]]-res_upwind_t[i].all_time[x_step_t,i1[j]])
		end
    end
	lable_t=["t=$ξ" for i in 1:1,ξ in time ]
	successive_Error=scatter(plot_upwind_t, ν_Range[1:end-1], err_t[:,2:end], label=lable_t )

	ploupwind_slope_t=plot(xlabel="ν",ylabel="Slope",title="slution convergence $method_name, x=$((x_step_t-1)*1/(nRange1[1]-1))*L ",legend= :bottomright)
    for i in 2:n_time
		plot!(ploupwind_slope_t,ν_Range[1:end-2], slope(err_t[:,i],ν_Range[1:end-1]),label="t=$(time[i])",lw=2)
    end
	return (slope_t=ploupwind_slope_t , successive_Error=successive_Error)
end
	

# ╔═╡ 790bfa26-e1d7-4ad8-91f0-3178bc84ba27
plot_ν_seris(time, ν_Range, res_upwind_t, x_step_t, nRange1,"upwind")

# ╔═╡ 11b1bf81-f368-4e69-a8c7-87c53dcdebec
"""
<div div style="direction: rtl;padding:14px">  
<p><b> شرایط اولیه نوع اول (n ثابت، ν متغیر ) </b>  </p>
<p>خطا برای شرایط اولیه نوع  اول :  </p>
</div>
"""|>HTML

# ╔═╡ dd17baaa-9dda-4f6e-99db-197fed04b820
res_upwind_t_err1= [ generate_all_time_upwind(n,inational_value!(L,n,1),L, a, νi,end_time,bc) for νi in ν_Range ]

# ╔═╡ 87030fac-e772-4bad-9827-56d9f7b1cd43
plot_ν_seris(time, ν_Range, res_upwind_t_err1, x_step_t, nRange1,"upwind")

# ╔═╡ a9c7c94c-2aa8-4c6b-b7d2-6dd099a27a24
"""
<div div style="direction: rtl;padding:14px">  
<p><b> شرایط اولیه نوع سوم (n ثابت، ν متغیر ) </b>  </p>
<p>خطا برای شرایط اولیه نوع  سوم:  </p>
</div>
"""|>HTML

# ╔═╡ 2a0a62a9-1640-42fb-ab92-67ffbfe4cafd
res_upwind_t_err3= [ generate_all_time_upwind(n,inational_value!(L,n,3),L, a, νi,end_time,bc) for νi in ν_Range ]

# ╔═╡ b1f12ce5-cc92-4b1e-8e4d-d37924b084be
plot_ν_seris(time, ν_Range, res_upwind_t_err3, x_step_t, nRange1,"upwind")

# ╔═╡ b7684934-c28f-4d12-bab5-7b521adfe9bb
"""<div style="direction: rtl"> <h2>   انجام محاسبات با روش Lax</h2><hr </div>

<div div style="direction: rtl;padding:14px">  
<p> از انجایی که مطالب توضیح داده شده در روش محاسبه FTCS بسیاری با سایر روش ها مشابه هستند تنها به ذکر مواتفاوت ها میپردازیم .  </p>
</div>


"""|>HTML

# ╔═╡ 3d15ed85-7eb1-4997-b5c6-79696b8d69be
"""<span style="direction: rtl"> <h3> تجزیه معادلات حاکم به روش  Lax</h3></span>

"""|>HTML

# ╔═╡ 581122a6-9190-4d0c-9ebe-875709f0c86c
md" 

```math
\begin{equation}
u_{i}^{ n+1}= \frac{(1+\nu)}{2}\times u_{i-1}^{n}+

\frac{(1-\nu)}{2}\times u_{i+1}^{n}
\end{equation}

```
"

# ╔═╡ e8538bd1-31b2-43ae-81a7-55a9cfbd3938


# ╔═╡ 8781a874-95bb-4352-a31c-bcd0ffa97845
md" 

```math
\begin{equation}
G_{m}= cos(\beta)+
i\nu\times sin(\beta)
\end{equation}

```
"

# ╔═╡ 60e52293-7d03-49d0-8604-a384d675f77b
begin
G_m_Lax(a,b)=cos(b)+(a*sin(b))im
#b=range(0, stop=π, length=100)
plo_G_m_Lax=scatter( title="polar plot G_m Lax",legend= :bottomright, proj=:polar, rlims=(0,1.5),reuse=true)
for a=1:-0.25:0.25
    #println(a)
    #θ=nothing
    #G_m=nothing
    θ=@.angle(G_m_Lax(a,b))
    G_m=@.abs(G_m_Lax(a,b))
    plot!(plo_G_m_Lax,θ, G_m,label="ν=$a")
end
	
plot(plo_G_m_Lax)
end

# ╔═╡ ff12be69-42da-455c-b948-a635592f0fdd
"""<span style="direction: rtl"> <h3> کد ها و مراحل محاسبات کامپیوتری روش Lax  </h3></span>

"""|>HTML

# ╔═╡ d1357278-e475-4e85-8491-306847558306
function internal_coeffs_Lax!(I, J, V, n, ν)
	for i= 2:(n-1)
		K=3*(i-2).+(1:3)# [3*n+1 ][3*o+1, 3*o+2, 3*o+3]
		I[K] .= i#[i, i, i]
		J[K] = [i-1, i, i+1]
		V[K] = [(1+ν)/2, 0, (1-ν)/2]
	end
end

# ╔═╡ 88b7f210-0324-43c0-8789-4a6038fd9bfc
function generate_les_Lax(n, ν, bc)
    h=dh(L,n)
	lbc,rbc=bc
	I, J, V= initilize_variables(n)
	internal_coeffs_Lax!(I, J, V, n, ν)
	lbc(I, J, V, n)
	rbc(I, J, V, n, ν)
	return sparse(I, J, V) #, C
end

# ╔═╡ 55b27c1d-8ae1-49b6-986a-eab21ea2383d
A_Lax = generate_les_Lax(n, ν, bc)

# ╔═╡ 71390fa2-e25b-4039-89ef-9f915e0b68d3
det(A_Lax)

# ╔═╡ 9c5d55aa-6c4a-40be-aff2-7b3bc3c9d347
function generate_all_time_Lax(n, u0, L, a, ν, end_time, bc)
	h_range=LinRange(0,L,n)
	τ=ν*h_range[2]/a
    time0=u0
	A = generate_les_Lax(n, ν, bc)
	time_lapse=collect(τ:τ:end_time)
	cloum=length(time_lapse)
	all_time=zeros(Float64,n,cloum)
	for i in (1:cloum)
		u1=A*(time0)
		all_time[:,i]=u1
		time0=u1
	end
	return (all_time=all_time , time_lapse=time_lapse, h_range=h_range)
end

# ╔═╡ 7880a698-2149-46bf-9010-b55bf32d22fb
u1_Lax=generate_all_time_Lax(n,inational_value!(L,n,1),L, a, ν,end_time,bc)

# ╔═╡ 280e11f6-4edb-449d-a873-91ed7161dabb
plo_lax=scatter(u1_Lax.h_range,u1_Lax.all_time[:,2],xlabel="ξ",ylabel="θ", label="n=$n")

# ╔═╡ 52370a1f-a76f-473c-8497-a7875e87790a


# ╔═╡ c69304fd-157f-4592-8454-38f3ca3716ae
"""
<div style="direction: rtl"> <h3>   تحلیل خطا در روش Lax</h3><hr </div>
<p>
<b>
ν ثابت، n متغیر
</b>
</p>

<div div style="direction: rtl;padding:14px"> <p> با توجه به نمودار Gm این روش ناپایداری است و انتظاری برای پایداری از آن نیست. این مورد در نتایج نیزمشهود است


   </p>



<p> </p>
</div>
"""|>HTML

# ╔═╡ b8a0a1a8-d7de-424f-b69f-938c441ce003
"""
<div div style="direction: rtl;padding:14px">  
<p><b> شرایط اولیه نوع دوم (ν ثابت، n متغیر) </b>  </p>
<p>خطا برای شرایط اولیه نوع  دوم:  </p>
</div>
"""|>HTML

# ╔═╡ 9b79090a-28b1-4f98-9513-cddd23b66c4f
res_Lax= [generate_all_time_Lax(ni,inational_value!(L,ni,2),L, a, ν,end_time,bc) for ni in nRange1 ]

# ╔═╡ 60cd1e9e-fed0-4cd0-bddc-5fd71d13914b
plot_x_seris_Lax=plot_x_seris(ν, res_Lax, time_want, nRange1,"Lax")

# ╔═╡ b634a9bd-9c01-4ddd-8f7b-ffea3305e12c
"""
<div div style="direction: rtl;padding:14px">  
<p> نمایش خطای نوع دوم برای n موجود در مجموعه nRange:</p>
</div>
"""|>HTML

# ╔═╡ 2e8b1518-1c61-47ee-b8e9-42dc3240a53d
plot_x_seris_Lax.solution

# ╔═╡ 4d4a6d93-8093-461c-a24d-bec5147f9391
plot_x_seris_Lax.successive_Error

# ╔═╡ 6771cdea-293d-4f68-b419-69511c89a322
plot_x_seris_Lax.plot_slope

# ╔═╡ a6992c0b-fa91-479c-90da-ad8b76bdc1c7
"""
<div div style="direction: rtl;padding:14px">  
<p><b> شرایط اولیه نوع اول (ν ثابت، n متغیر)  </b>  </p>
<p>خطا برای شرایط اولیه نوع  اول :  </p>
</div>
"""|>HTML

# ╔═╡ 5677de1a-ada9-4924-9d05-4a24c7144b59
res_Lax_err_1= [generate_all_time_Lax(ni,inational_value!(L,ni,1),L, a, ν,end_time,bc) for ni in nRange1 ]

# ╔═╡ 7b007d76-8485-41d0-98c5-a7e03383797b
plot_x_seris(ν, res_Lax_err_1, time_want, nRange1,"Lax")

# ╔═╡ f57f9d5f-a2b7-47b8-bd20-e1b6ddee3e2f
"""
<div div style="direction: rtl;padding:14px">  
<p><b> شرایط اولیه نوع سوم (ν ثابت، n متغیر)  </b>  </p>
<p>خطا برای شرایط اولیه نوع  سوم:  </p>
</div>
"""|>HTML

# ╔═╡ 72312a27-f66d-40a5-8a02-214d7e0115e9
res_Lax_err_3= [generate_all_time_Lax(ni,inational_value!(L,ni,3),L, a, ν,end_time,bc) for ni in nRange1 ]

# ╔═╡ 203b145c-f255-4775-b8b1-8a095b5bd933
plot_x_seris(ν, res_Lax_err_3, time_want, nRange1,"Lax")

# ╔═╡ 8e0b6052-1214-4d1f-afe6-69183c07056e
"""
<div div style="direction: rtl;padding:14px">  
<p> 
<b>
n ثابت، ν متغیر 
</b>
</p>
<p> 
در این بخش با ثابت نگه داشتن ν و تغییر n در واقع τ را تغییر میدهیم.این نمودار های با معنی تر از نمودار های قبلی هستند و بهتر تفصیر میشوند.
</p>
</div>
"""|>HTML

# ╔═╡ e25349c9-646c-45f3-a973-6a93a8530891
"""
<div div style="direction: rtl;padding:14px">  
<p><b> شرایط اولیه نوع دوم (n ثابت، ν متغیر ) </b>  </p>
<p>خطا برای شرایط اولیه نوع  دوم:  </p>
</div>
"""|>HTML

# ╔═╡ 570bd5fa-dea4-402c-af2b-2ded5a56bb26
res_Lax_t= [ generate_all_time_Lax(n,inational_value!(L,n,2),L, a, νi,end_time,bc) for νi in ν_Range ]

# ╔═╡ ac2eca55-954c-4d5b-a1dd-f2eb9be341e7
plot_ν_seris_Lax=plot_ν_seris(time, ν_Range, res_Lax_t, x_step_t, nRange1,"Lax")

# ╔═╡ 074264b8-4781-4d0a-b406-14278484332a
plot_ν_seris_Lax.slope_t

# ╔═╡ 3dc7ffe4-4eec-4af0-b175-5489e1c89b27
plot_ν_seris_Lax.successive_Error

# ╔═╡ 1189492c-d978-4a23-9241-42a99d8fa1cd
"""
<div div style="direction: rtl;padding:14px">  
<p><b> شرایط اولیه نوع اول (n ثابت، ν متغیر ) </b>  </p>
<p>خطا برای شرایط اولیه نوع  اول :  </p>
</div>
"""|>HTML

# ╔═╡ ee8dbdb6-854b-4efb-b555-4dca2c69b035
res_Lax_t_err1= [ generate_all_time_Lax(n,inational_value!(L,n,1),L, a, νi,end_time,bc) for νi in ν_Range ]

# ╔═╡ 6555c135-b218-4e4c-9b88-c0e40cd692d4
plot_ν_seris(time, ν_Range, res_Lax_t_err1, x_step_t, nRange1,"Lax")

# ╔═╡ 2420ddad-f9c7-4328-bae1-780939c11287
"""
<div div style="direction: rtl;padding:14px">  
<p><b> شرایط اولیه نوع سوم (n ثابت، ν متغیر ) </b>  </p>
<p>خطا برای شرایط اولیه نوع  سوم:  </p>
</div>
"""|>HTML

# ╔═╡ 9aa80b00-6b21-4af3-80d5-6f28ae9a46fe
res_Lax_t_err3= [ generate_all_time_Lax(n,inational_value!(L,n,3),L, a, νi,end_time,bc) for νi in ν_Range ]

# ╔═╡ e239f9ef-1a13-4cd9-bc45-1d677d57e2fb
plot_ν_seris(time, ν_Range, res_Lax_t_err3, x_step_t, nRange1,"Lax")

# ╔═╡ a030b0f5-6d4a-415c-99f5-dc65c56ef90b
"""<div style="direction: rtl"> <h2>   انجام محاسبات با روش: Lax-Wendroff</h2> </div>
<div div style="direction: rtl;padding:14px">  
<p> روش CTCS به دلیل تفاوت در تجریه و معماری با سایر روش ها تفاوت جزئی دارد(نیاز به دو گام زمانی قبل برای محاسبه گام بعدی) و این موضوع باعث شده تا نیاز به باز تعریف تعدادی از توابقع اولیه داشته باشیم (البته زیاد لزوم نداشت اما از این که توابع قبلی را جامع تر کنیم و تغییر سراسری بدیم راحت تر بود )   </p>
</div>

"""|>HTML

# ╔═╡ 1ef6f8b0-6e6f-4b0b-976d-1a9ae8865768
"""<span style="direction: rtl"> <h3> تجزیه معادلات حاکم به روش  Lax-Wendroff</h3></span>

"""|>HTML

# ╔═╡ 593c689f-2277-46c2-88bc-8949e71b7965
md" 

```math
\begin{equation}
u_{i}^{ n+1}= \frac{(\nu^2+\nu)}{2}\times u_{i-1}^{n}+
(1-\nu^2)\times u_{i}^{ n}+
\frac{(\nu^2-\nu)}{2}\times u_{i+1}^{n}
\end{equation}

```
"

# ╔═╡ 7b0e6aea-a56e-4a48-8941-9e8d4cd97e19
md" 

```math
\begin{equation}
G_{m}= 1-\nu^2+\frac{\nu^2}{2}\times cos(\beta)+
i\nu\times sin(\beta)
\end{equation}

```
"

# ╔═╡ 694ce6ca-2e0a-4b87-8d41-e44275a15306
"""
<div div style="direction: rtl;padding:14px">  



<p> 
با توجه به نمودار Gm این روش تنها برای ν های کوچک پایداری است. این مورد در نتایج نیزمشهود است. با کوچک شدن بیشتر ν پایدار تر میشود. همچمین اختلاف فاز هارمونیک ها نیز کاهش میابد  .
</p>

</div>
"""|>HTML



# ╔═╡ fe2822b1-d2be-484c-9b27-853fa1aa9139
begin
G_m_Lax_Wendroff(a,b)=1-a^2+(a^2/2)*cos(b)+(a*sin(b))im
#b=range(0, stop=π, length=100)
plo_G_m_Lax_Wendroff=plot(title="polar plot G_m Lax-Wendroff",legend= :bottomright, proj=:polar, rlims=(0,1.5),reuse=true)
for a=1:-0.25:0.25
    #println(a)
    #θ=nothing
    #G_m=nothing
    θ=@.angle(G_m_Lax_Wendroff(a,b))
    G_m=@.abs(G_m_Lax_Wendroff(a,b))
    plot!(plo_G_m_Lax_Wendroff,θ, G_m,label="ν=$a")
end
	
plot(plo_G_m_Lax_Wendroff)
end

# ╔═╡ d3bdf9e9-88e8-46e8-95de-938634e6e7a3
"""<span style="direction: rtl"> <h3> کد ها و مراحل محاسبات کامپیوتری روش Lax-Wendroff  </h3></span>

"""|>HTML

# ╔═╡ 1121bdcb-3dda-49ad-a05d-7c30c3a05f65
function internal_coeffs_Lax_Wendroff!(I, J, V, n, ν)
	for i= 2:(n-1)
		K=3*(i-2).+(1:3)# [3*n+1 ][3*o+1, 3*o+2, 3*o+3]
		I[K] .= i#[i, i, i]
		J[K] = [i-1, i, i+1]
		V[K] = [(ν^2+ν)/2, (1-ν^2), (ν^2-ν)/2]
	end
end

# ╔═╡ 93f2d1c9-ad73-49f5-8fd1-f4bda61278fa
function generate_les_Lax_Wendroff(n, ν, bc)
    h=dh(L,n)
	lbc,rbc=bc
	I, J, V= initilize_variables(n)
	internal_coeffs_Lax_Wendroff!(I, J, V, n, ν)
	lbc(I, J, V, n)
	rbc(I, J, V, n, ν)
	return sparse(I, J, V) #, C
end

# ╔═╡ 31837f14-7acc-4b4d-bb17-feecb9bda924
A_Lax_Wendroff = generate_les_Lax_Wendroff(n, ν, bc)

# ╔═╡ 01224627-3227-4c82-aba7-78e8e8a8610c
function generate_all_time_Lax_Wendroff(n, u0, L, a, ν, end_time, bc)
	h_range=LinRange(0,L,n)
	τ=ν*h_range[2]/a
    time0=u0
	A = generate_les_Lax_Wendroff(n, ν, bc)
	time_lapse=collect(τ:τ:end_time)
	cloum=length(time_lapse)
	all_time=zeros(Float64,n,cloum)
	for i in (1:cloum)
		u1=A*(time0)
		all_time[:,i]=u1
		time0=u1
	end
	return (all_time=all_time , time_lapse=time_lapse, h_range=h_range)
end

# ╔═╡ eb38130c-a2f9-4c4e-b30d-19cb12880e0e
u1_Lax_Wendroff=generate_all_time_Lax_Wendroff(n,inational_value!(L,n,1),L, a, ν,end_time,bc)

# ╔═╡ 7ecebffa-9a96-43a0-928d-d3b943879452
"""
<div style="direction: rtl"> <h3>   تحلیل خطا در روش Lax_Wendroff</h3> </div>



<div div style="direction: rtl;padding:14px"> 
<p>
<b>
ν ثابت، n متغیر
</b>
</p>

<div div style="direction: rtl;padding:14px"> <p> با توجه به نمودار Gm این روش تنها برای ν های کوچک پایداری است. این مورد در نتایج نیزمشهود است.


   </p>


<p> </p>
</div>
"""|>HTML

# ╔═╡ e1556781-6273-4277-b8ae-d6b1adea8508
"""
<div div style="direction: rtl;padding:14px">  
<p><b> شرایط اولیه نوع دوم (ν ثابت، n متغیر) </b>  </p>
<p>خطا برای شرایط اولیه نوع  دوم:  </p>
</div>
"""|>HTML

# ╔═╡ 1f5246a6-e9db-4319-a9a2-ee0182007d43
 res_Lax_Wendroff= [generate_all_time_Lax_Wendroff(ni,inational_value!(L,ni,2),L, a, 0.5*ν,end_time,bc) for ni in nRange1 ]

# ╔═╡ ea73204c-42b8-43a7-a2b6-5e82fae660e1
plot_x_seris_Lax_Wendroff=plot_x_seris(0.5*ν, res_Lax_Wendroff, time_want, nRange1,"Lax_Wendroff")

# ╔═╡ ceb4fd5e-805f-41c4-a9d8-299dde6c4331
plot_x_seris_Lax_Wendroff.solution

# ╔═╡ b8c4cc42-81b6-46fe-a947-b211215086d4
plot_x_seris_Lax_Wendroff.successive_Error

# ╔═╡ 75ee1183-05d1-4644-93a0-ac007d909c61
"""
<div div style="direction: rtl;padding:14px">  

<p>همان طور که قابل مشاهده است ثابت ماندن ν تاثیر زیادی در کاهش استهلاک ندارد و در نواحی مرزی مثل مرز سمت چپ ناپایداری مشاهده میشود(نقاط x<0.4) اما به مرور هر چه به سمت مرکز شبکه)(x>0.4) میرویم همگرایی  با کاهش h قابل مشاهده است.  </p>
</div>
"""|>HTML

# ╔═╡ de6c30a0-d8be-4772-a3df-ad1dc5737e9c
plot_x_seris_Lax_Wendroff.plot_slope

# ╔═╡ 396755b4-3d27-4d67-a519-9b40641f5619
"""
<div div style="direction: rtl;padding:14px">  
<p><b> شرایط اولیه نوع اول (ν ثابت، n متغیر)  </b>  </p>
<p>خطا برای شرایط اولیه نوع  اول :  </p>
</div>
"""|>HTML

# ╔═╡ c376fb47-b1ac-4383-b501-325a5075778b
res_Lax_Wendroff_err_1= [generate_all_time_Lax_Wendroff(ni,inational_value!(L,ni,1),L, a, ν,end_time,bc) for ni in nRange1 ]

# ╔═╡ 6843544d-4678-4328-877b-2e842f2ddb65
plot_x_seris(ν, res_Lax_Wendroff_err_1, time_want, nRange1,"Lax-Wendroff")

# ╔═╡ 9d3a6f9f-049f-4440-b0dd-3424f8a24abb
"""
<div div style="direction: rtl;padding:14px">  
<p><b> شرایط اولیه نوع سوم (ν ثابت، n متغیر)  </b>  </p>
<p>خطا برای شرایط اولیه نوع  سوم:  </p>
</div>
"""|>HTML

# ╔═╡ aaf9ed7d-d06d-4d93-913f-9c228ba80b2d
res_Lax_Wendroff_err_3= [generate_all_time_Lax_Wendroff(ni,inational_value!(L,ni,3),L, a, ν,end_time,bc) for ni in nRange1 ]

# ╔═╡ 76c3bc04-2357-47f9-a1f3-14f82a05ee04
plot_x_seris(ν, res_Lax_Wendroff_err_3, time_want, nRange1,"Lax-Wendroff")

# ╔═╡ af3c2185-9948-45eb-9b23-918d8002c769
"""
<div div style="direction: rtl;padding:14px">  
<p> 
<b>
n ثابت، ν متغیر 
</b>
</p>
<p> 
در این بخش با ثابت نگه داشتن ν و تغییر n در واقع τ را تغییر میدهیم.این نمودار های با معنی تر از نمودار های قبلی هستند و بهتر تفصیر میشوند.
</p>
</div>
"""|>HTML

# ╔═╡ 842ab124-7517-433d-bf11-bf88cb583b2d
"""
<div div style="direction: rtl;padding:14px">  
<p><b> شرایط اولیه نوع دوم (n ثابت، ν متغیر ) </b>  </p>
<p>خطا برای شرایط اولیه نوع  دوم:  </p>
</div>
"""|>HTML

# ╔═╡ ed8e432e-5080-4c77-bce7-9c6b1117e437
res_Lax_Wendroff_t= [ generate_all_time_Lax_Wendroff(n,inational_value!(L,n,2),L, a, νi,end_time,bc) for νi in ν_Range ]

# ╔═╡ 198b4615-8e34-4d8f-808f-7df7f6e7bd35
plot_ν_seris_Lax_Wendroff=plot_ν_seris(time, ν_Range, res_Lax_Wendroff_t, x_step_t, nRange1,"Lax")

# ╔═╡ 0a3b8c46-5c5d-4f81-8df0-4900ee8cc140
plot_ν_seris_Lax_Wendroff.successive_Error

# ╔═╡ 1588e2fc-1716-40da-bba4-985c693e9176
plot_ν_seris_Lax_Wendroff.slope_t

# ╔═╡ 208a37b3-e886-42f4-9302-53968d7ed446
"""
<div div style="direction: rtl;padding:14px">  
<p><b> شرایط اولیه نوع اول (n ثابت، ν متغیر ) </b>  </p>
<p>خطا برای شرایط اولیه نوع  اول :  </p>
</div>
"""|>HTML

# ╔═╡ e04cba86-6c99-4f50-857b-4f530f1efee2
res_Lax_Wendroff_t_err1= [ generate_all_time_Lax_Wendroff(n,inational_value!(L,n,1),L, a, νi,end_time,bc) for νi in ν_Range ]

# ╔═╡ 63f9e66d-799e-4d6e-9833-d9b12275b473
plot_ν_seris(time, ν_Range, res_Lax_Wendroff_t_err1, x_step_t, nRange1,"Lax-Wendroff")

# ╔═╡ bebdc082-b5f1-47fe-bf1b-e9dd93306d47
"""
<div div style="direction: rtl;padding:14px">  
<p><b> شرایط اولیه نوع سوم (n ثابت، ν متغیر ) </b>  </p>
<p>خطا برای شرایط اولیه نوع  سوم:  </p>
</div>
"""|>HTML

# ╔═╡ a063c3fc-d716-44e3-ab1a-9af708b915a0
res_Lax_Wendroff_t_err3= [ generate_all_time_Lax_Wendroff(n,inational_value!(L,n,3),L, a, νi,end_time,bc) for νi in ν_Range ]

# ╔═╡ bdca1cd9-b5d1-4e7f-9c61-73f344a55988
plot_ν_seris(time, ν_Range, res_Lax_Wendroff_t_err3, x_step_t, nRange1,"Lax-Wendroff")

# ╔═╡ 880bc247-7440-41f0-89d1-28ea6874fae8
"""<div style="direction: rtl"> <h2>   انجام محاسبات با روش: MacCormack</h2><hr </div>


"""|>HTML

# ╔═╡ 6751ebaa-56a5-4312-96cb-31fe084ef7ae
"""<span style="direction: rtl"> <h3> تجزیه معادلات حاکم به روش MacCormack</h3></span>

"""|>HTML

# ╔═╡ fe70925e-c122-44fe-965f-74808407c617
"""
<div div style="direction: rtl;padding:14px">  
<p><b> کاملا مشابه روش Lax-Wendroff </b>  </p>
<p>  </p>
</div>
"""|>HTML

# ╔═╡ 33fe3610-e5bf-4169-9bb8-8af757d591fa
md" 

```math
\begin{equation}
u_{i}^{ n+1}= \frac{(\nu^2+\nu)}{2}\times u_{i-1}^{n}+
(1-\nu^2)\times u_{i}^{ n}+
\frac{(\nu^2-\nu)}{2}\times u_{i+1}^{n}
\end{equation}

```
"

# ╔═╡ aa6a7f30-72e3-4ec7-b325-f470fa2d2682
md" 

```math
\begin{equation}
G_{m}= 1-\nu^2+\frac{\nu^2}{2}\times cos(\beta)+
i\nu\times sin(\beta)
\end{equation}

```
"

# ╔═╡ 7bd41cf4-ed3c-40fd-8ac1-86f4bb448f22
"""<span style="direction: rtl"> <h3> کد ها و مراحل محاسبات کامپیوتری روش MacCormack </h3>
<p> این روش در معادلات خطی (مانند معادله حل شده ما) پیکره محاسباتی کاملا مشابه روش Lax-Wendroff دارد(همان طور که در تجزیه قابل مشاهده است) و در محاسبات بین این دو روش تفاوتی نیست و از تکرار آن خودداری شده است .   </p>
</span>

"""|>HTML

# ╔═╡ 84de0ab9-d481-47b9-9229-bbeebae3b45a
"""
<div div style="direction: rtl"> 
<p>
</p>
<p> </p>
</div>
"""|>HTML

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[compat]
Plots = "~1.27.1"
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
git-tree-sha1 = "91b5dcf362c5add98049e6c29ee756910b03051d"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.3"

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
git-tree-sha1 = "4f00cc36fede3c04b8acf9b2e2763decfdcecfa6"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.13"

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
git-tree-sha1 = "58f25e56b706f95125dcb796f39e1fb01d913a71"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.10"

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
git-tree-sha1 = "ab05aa4cc89736e95915b01e7279e61b1bfe33b8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.14+0"

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
git-tree-sha1 = "bb16469fd5224100e422f0b027d26c5a25de1200"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.2.0"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "1690b713c3b460c955a2957cd7487b1b725878a7"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.27.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "d3538e7f8a790dc8903519090857ef8e1283eecd"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.5"

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
git-tree-sha1 = "6976fab022fea2ffea3d945159317556e5dad87c"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.2"

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
# ╟─35e7ba37-23a2-4d03-8be0-b01997c71c55
# ╟─e21545ec-f882-475e-ac99-2e4b3d6525c5
# ╟─0c71ef7b-dbaa-4918-a2ad-99cc49500757
# ╟─fe026e1a-8bee-46cb-8f0a-f9e390321f23
# ╠═56eee02f-091a-4b5c-a850-08fb59f12aae
# ╠═9a1a0057-54a5-4c33-93ff-c3cdd70f1662
# ╟─aece47b6-9aa9-474d-ab82-da97eb53c2ad
# ╠═ee228b61-0538-49f0-ab30-74246458cf51
# ╟─5be587ad-7703-4430-a9eb-279a24d20fc0
# ╠═68b3dff3-1a99-4c31-ad3d-3b2455de69a6
# ╠═c34595ea-f229-47e8-aa01-b187bc9f90c2
# ╠═b90a98a0-e085-4f36-b415-c89a0bff1624
# ╠═5a4bc4e6-1bad-4a4a-bd73-00d56055d895
# ╠═94915014-48b1-4c3c-93bb-5f467c3d3502
# ╟─28fb4806-a658-4489-9296-c36b1e1c504e
# ╠═118ad3b0-fec2-4b86-9e58-c1f099602819
# ╟─60e0b122-25da-42a5-80e9-bf42bb6f6455
# ╠═39d600fd-a18d-4340-96af-2b2ccd3f50fc
# ╟─e4c74479-e565-4570-bd3a-fc4ec8c52b46
# ╟─c2df185f-5284-41e9-a662-8246b27d684b
# ╟─509b28b6-5152-4dc3-b2b8-4ec4a2101df8
# ╟─9819c30a-08c8-4ca4-9f4b-e158b35a89d8
# ╟─cfbc4bb5-56bb-4659-9036-71e7e092744f
# ╠═14146054-6afb-456c-8db9-94fcde00d331
# ╟─1c68b2a3-62c6-4449-9e2c-6c69b9a3df04
# ╠═ec51d580-7305-4be4-9e64-6507c1230cac
# ╠═eb548fea-9c56-4b39-9d2a-2caa9926b80f
# ╠═dd23acef-14fe-41c3-b991-b61256b24d75
# ╠═1c3a2a15-a995-4124-9bf2-c61b5c62af23
# ╟─4453c8b4-e94c-4095-a3cf-4ac64be7d35d
# ╟─bfc059e1-2561-46d9-b2c1-71adbdbdbc7e
# ╠═022e31c3-8864-4106-9aa4-6ebd49124814
# ╟─b5556da7-056d-493f-bded-40d73376bd58
# ╟─b3efc4cc-25a0-46a7-b1e4-7d6d15194f32
# ╟─ae7ba8c0-d12c-45c2-99f4-863b52f314b7
# ╟─1017451b-46c3-4009-8fd8-ba257db39ac2
# ╠═a2daaecd-b154-4363-9feb-2c87d3300538
# ╠═65de7f54-9384-4124-b353-fa1ef495a1a6
# ╟─9c060b1d-a057-489d-844f-69992df65bc1
# ╟─fef9cd8b-05a2-403e-a092-81d9960feeff
# ╟─503388cf-eeb1-4c42-b04f-0263bcfb316e
# ╠═9e25af28-7ba6-4db1-a659-36e5e43aeea1
# ╠═287d567a-7905-4801-a671-be58fcd5b8b0
# ╠═051ec1b2-33fe-4119-925a-ec1a34a20f63
# ╟─2201cf8f-8e48-47db-a5f0-f6cc2265ef8f
# ╠═3b7e65df-07de-4c5c-8e42-e7159aeeb752
# ╠═1f3db9fb-0156-4f7d-bd6c-14fd4f512a45
# ╟─866d2f41-fdc3-4e49-ba91-49d95999eeed
# ╟─9d190dbf-afcf-4293-97c5-13d2d5b7e5e8
# ╟─c76ae0c1-ed35-4547-9fef-8f3c8a464092
# ╟─6743e0ca-b724-4822-94df-e8254853c842
# ╟─9b33684a-337d-424b-a632-881fbe4e7571
# ╟─007f9a43-f4c8-4d3e-96f9-895de84efa97
# ╠═a56a5e2f-4302-438a-8921-c31e13674b22
# ╟─6d92754d-c4f0-4b34-ad0b-9f1b57391331
# ╟─8fc18020-e7bf-4060-abcb-5a4ff3c636df
# ╟─76164d36-37e4-4d6b-8d7b-b315c1c5c801
# ╠═be2591de-b5f1-4e16-b7d9-56c864c362db
# ╟─7f1cd37b-27fb-45fb-961e-9bfd2e979b55
# ╠═fc63e35c-0b20-4a36-8760-543db4f0cac5
# ╠═a68b1cca-7642-4593-8ba0-04b36931c20a
# ╠═70ea54df-7d89-4411-90a4-56230db9c868
# ╟─3a950fd4-b631-4868-bcb7-07089d2eb4f9
# ╠═0b754307-067b-4ab8-b857-7b340c82eb50
# ╠═e2649731-93ce-4783-942f-82b6954bf061
# ╠═b2398175-5d89-4311-9651-0b494984a37a
# ╟─befb0951-9cc9-4f5c-ac6d-81ae08842b66
# ╠═fdc80205-b52b-4606-a6f9-16288c674299
# ╠═1878e5d0-ff15-4a6d-a2da-30f86bf9ab22
# ╟─2b7a7b49-97f8-4c02-aa0f-b855bc90b268
# ╟─e218cf86-2e2b-412f-b6ed-e5d6ba340f39
# ╠═f68d5d7e-3f64-4b04-9aab-31448a59fd3c
# ╟─f069dfeb-0cc6-49df-86c3-5aa60f807f99
# ╠═583a8f2c-913b-4eda-ad90-f66b3a38d87b
# ╟─3689c5df-404d-483d-91a1-c3df537deb62
# ╠═0a9f70c5-9be2-44d3-8541-7f0111977ca6
# ╠═2d92c93f-4e09-48a3-b92a-fb6bc9d49d24
# ╟─561a0b16-6865-4381-bb48-9dfa79e35c37
# ╠═ec900552-34c8-4486-b37c-a32536f61cfb
# ╟─6f2f21c1-e9a7-4db4-a1c2-6204abe2ae48
# ╠═ef5ddb37-91fd-446b-8333-e50a48008fc2
# ╠═e5659b9e-5d9a-43ed-8a67-3b2245ffb138
# ╠═52e072cd-a21e-4a83-88b8-3fe3678c0fd8
# ╟─0e5f9639-c5a6-4a17-b077-8ad37361de80
# ╠═ba0a7ab7-c687-4a87-9619-2e3e1730199a
# ╠═467bf194-1daa-44d4-a1c7-78863d293555
# ╠═ef11bf7a-33b6-4ef1-afa6-3077e5efd9ac
# ╟─47ac2a22-20c7-4a43-95ae-f9f7195a71a1
# ╠═07cb20d4-1571-4bec-9dae-00017e03ce07
# ╠═5e759c6c-2ce9-4c1d-89de-01a700b015de
# ╟─94dce462-9662-402f-85d3-074df56d96ab
# ╠═8c85bfb2-d5bd-45d5-8368-ab25f94ebc45
# ╠═4699894f-4560-4658-9b4a-960eb7fb6b5f
# ╠═8a790d36-de2f-485c-b476-fee78772c796
# ╠═3d874be8-61c8-4620-bd2b-fbdb00a9f3ab
# ╠═2bb874d6-2018-4c07-930f-d040bcc5eba7
# ╠═ae7a2f50-ae54-49fa-8119-0ae062ee1a30
# ╠═dc43d2af-2975-4c94-80aa-c0cf7c189554
# ╠═201ff5bc-8d24-4400-bcbc-d16262c22347
# ╟─061d2b74-96e1-43b8-b252-13911ba9639a
# ╠═89aea175-f4e1-4777-ba28-0c684a959295
# ╟─27883575-5235-40b1-97aa-d1779b5f35b6
# ╠═76fac5ab-9676-4bfe-a08b-4ce75ea98774
# ╠═790bfa26-e1d7-4ad8-91f0-3178bc84ba27
# ╟─11b1bf81-f368-4e69-a8c7-87c53dcdebec
# ╠═dd17baaa-9dda-4f6e-99db-197fed04b820
# ╠═87030fac-e772-4bad-9827-56d9f7b1cd43
# ╟─a9c7c94c-2aa8-4c6b-b7d2-6dd099a27a24
# ╠═2a0a62a9-1640-42fb-ab92-67ffbfe4cafd
# ╠═b1f12ce5-cc92-4b1e-8e4d-d37924b084be
# ╟─b7684934-c28f-4d12-bab5-7b521adfe9bb
# ╟─3d15ed85-7eb1-4997-b5c6-79696b8d69be
# ╟─581122a6-9190-4d0c-9ebe-875709f0c86c
# ╟─e8538bd1-31b2-43ae-81a7-55a9cfbd3938
# ╟─8781a874-95bb-4352-a31c-bcd0ffa97845
# ╠═60e52293-7d03-49d0-8604-a384d675f77b
# ╟─ff12be69-42da-455c-b948-a635592f0fdd
# ╠═d1357278-e475-4e85-8491-306847558306
# ╠═88b7f210-0324-43c0-8789-4a6038fd9bfc
# ╠═55b27c1d-8ae1-49b6-986a-eab21ea2383d
# ╠═71390fa2-e25b-4039-89ef-9f915e0b68d3
# ╠═9c5d55aa-6c4a-40be-aff2-7b3bc3c9d347
# ╠═7880a698-2149-46bf-9010-b55bf32d22fb
# ╠═280e11f6-4edb-449d-a873-91ed7161dabb
# ╟─52370a1f-a76f-473c-8497-a7875e87790a
# ╟─c69304fd-157f-4592-8454-38f3ca3716ae
# ╠═b8a0a1a8-d7de-424f-b69f-938c441ce003
# ╠═9b79090a-28b1-4f98-9513-cddd23b66c4f
# ╠═60cd1e9e-fed0-4cd0-bddc-5fd71d13914b
# ╟─b634a9bd-9c01-4ddd-8f7b-ffea3305e12c
# ╠═2e8b1518-1c61-47ee-b8e9-42dc3240a53d
# ╠═4d4a6d93-8093-461c-a24d-bec5147f9391
# ╠═6771cdea-293d-4f68-b419-69511c89a322
# ╠═a6992c0b-fa91-479c-90da-ad8b76bdc1c7
# ╠═5677de1a-ada9-4924-9d05-4a24c7144b59
# ╠═7b007d76-8485-41d0-98c5-a7e03383797b
# ╠═f57f9d5f-a2b7-47b8-bd20-e1b6ddee3e2f
# ╠═72312a27-f66d-40a5-8a02-214d7e0115e9
# ╠═203b145c-f255-4775-b8b1-8a095b5bd933
# ╟─8e0b6052-1214-4d1f-afe6-69183c07056e
# ╠═e25349c9-646c-45f3-a973-6a93a8530891
# ╠═570bd5fa-dea4-402c-af2b-2ded5a56bb26
# ╠═ac2eca55-954c-4d5b-a1dd-f2eb9be341e7
# ╠═074264b8-4781-4d0a-b406-14278484332a
# ╠═3dc7ffe4-4eec-4af0-b175-5489e1c89b27
# ╟─1189492c-d978-4a23-9241-42a99d8fa1cd
# ╠═ee8dbdb6-854b-4efb-b555-4dca2c69b035
# ╠═6555c135-b218-4e4c-9b88-c0e40cd692d4
# ╟─2420ddad-f9c7-4328-bae1-780939c11287
# ╠═9aa80b00-6b21-4af3-80d5-6f28ae9a46fe
# ╠═e239f9ef-1a13-4cd9-bc45-1d677d57e2fb
# ╠═a030b0f5-6d4a-415c-99f5-dc65c56ef90b
# ╟─1ef6f8b0-6e6f-4b0b-976d-1a9ae8865768
# ╟─593c689f-2277-46c2-88bc-8949e71b7965
# ╟─7b0e6aea-a56e-4a48-8941-9e8d4cd97e19
# ╟─694ce6ca-2e0a-4b87-8d41-e44275a15306
# ╠═fe2822b1-d2be-484c-9b27-853fa1aa9139
# ╟─d3bdf9e9-88e8-46e8-95de-938634e6e7a3
# ╠═1121bdcb-3dda-49ad-a05d-7c30c3a05f65
# ╠═93f2d1c9-ad73-49f5-8fd1-f4bda61278fa
# ╠═31837f14-7acc-4b4d-bb17-feecb9bda924
# ╠═01224627-3227-4c82-aba7-78e8e8a8610c
# ╠═eb38130c-a2f9-4c4e-b30d-19cb12880e0e
# ╟─7ecebffa-9a96-43a0-928d-d3b943879452
# ╟─e1556781-6273-4277-b8ae-d6b1adea8508
# ╠═1f5246a6-e9db-4319-a9a2-ee0182007d43
# ╠═ea73204c-42b8-43a7-a2b6-5e82fae660e1
# ╠═ceb4fd5e-805f-41c4-a9d8-299dde6c4331
# ╠═b8c4cc42-81b6-46fe-a947-b211215086d4
# ╟─75ee1183-05d1-4644-93a0-ac007d909c61
# ╠═de6c30a0-d8be-4772-a3df-ad1dc5737e9c
# ╠═396755b4-3d27-4d67-a519-9b40641f5619
# ╠═c376fb47-b1ac-4383-b501-325a5075778b
# ╠═6843544d-4678-4328-877b-2e842f2ddb65
# ╟─9d3a6f9f-049f-4440-b0dd-3424f8a24abb
# ╠═aaf9ed7d-d06d-4d93-913f-9c228ba80b2d
# ╠═76c3bc04-2357-47f9-a1f3-14f82a05ee04
# ╟─af3c2185-9948-45eb-9b23-918d8002c769
# ╟─842ab124-7517-433d-bf11-bf88cb583b2d
# ╠═ed8e432e-5080-4c77-bce7-9c6b1117e437
# ╠═198b4615-8e34-4d8f-808f-7df7f6e7bd35
# ╠═0a3b8c46-5c5d-4f81-8df0-4900ee8cc140
# ╠═1588e2fc-1716-40da-bba4-985c693e9176
# ╠═208a37b3-e886-42f4-9302-53968d7ed446
# ╠═e04cba86-6c99-4f50-857b-4f530f1efee2
# ╠═63f9e66d-799e-4d6e-9833-d9b12275b473
# ╟─bebdc082-b5f1-47fe-bf1b-e9dd93306d47
# ╠═a063c3fc-d716-44e3-ab1a-9af708b915a0
# ╠═bdca1cd9-b5d1-4e7f-9c61-73f344a55988
# ╟─880bc247-7440-41f0-89d1-28ea6874fae8
# ╟─6751ebaa-56a5-4312-96cb-31fe084ef7ae
# ╟─fe70925e-c122-44fe-965f-74808407c617
# ╟─33fe3610-e5bf-4169-9bb8-8af757d591fa
# ╟─aa6a7f30-72e3-4ec7-b325-f470fa2d2682
# ╟─7bd41cf4-ed3c-40fd-8ac1-86f4bb448f22
# ╟─84de0ab9-d481-47b9-9229-bbeebae3b45a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
