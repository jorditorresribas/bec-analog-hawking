# Simulación de radiación de Hawking análoga con revisión cruzada entre dos IAs

**Un agujero negro acústico, un resultado espurio y cómo se identificó su origen**

**Autor:** Jordi Torres Ribas · **Contacto:** jorditorresribas+bec@gmail.com
**Versión:** 1.5 · **Fecha:** 18 de noviembre de 2026

---

El objetivo era simular en un ordenador un análogo de laboratorio de un agujero negro: un gas de átomos ultrafríos que fluye más rápido que el sonido y que, según la teoría, debería emitir una débil radiación térmica, el equivalente acústico de la radiación de Hawking. A mitad del trabajo, la simulación produjo un resultado llamativo que podría haberse presentado como hallazgo. No lo era. El documento incluye los errores cometidos por el camino, porque son la parte más informativa.

## Resumen

- Construí con IA un simulador de un análogo de laboratorio de un agujero negro y lo sometí a reglas de validación fijadas antes de cada prueba.
- A mitad del trabajo apareció un resultado llamativo: la radiación "colapsaba" al cambiar un parámetro del modelo. **No era física: era un error del modelo.** El test de simetría identificó el defecto: el residuo medido coincidía con el predicho en cuatro cifras significativas. Al corregirlo, el "colapso" desapareció.
- El instrumento reproduce un resultado **exacto** publicado por otro grupo de investigación: primero con un error de 3·10⁻⁵ en el coeficiente principal y, en una segunda prueba con los criterios publicados antes de ejecutarla, con un error de 10⁻⁸.
- Corregido el error, resultó que el parámetro que se estaba variando no hacía lo que se pretendía: equivale exactamente a cambiar la anchura del horizonte.
- Una segunda fase, para aislar de verdad la física de escalas pequeñas, terminó **sin veredicto**: un control fijado de antemano detectó que el observable elegido no medía lo que se pretendía, y detuvo una conclusión que los datos no justificaban.
- Las auditorías detectaron errores reales en los propios controles y en la documentación. Están documentados más abajo, con quién los detectó.

## Qué demuestra este proyecto

1. **Dirigir IA en un problema técnico difícil sin que los errores lleguen a las conclusiones**: una IA construye, otra audita, y una persona decide.
2. **Convertir "parece que funciona" en "está comprobado"**: criterios de éxito y de fracaso escritos *antes* de ejecutar. Hay dos excepciones, las dos en tolerancias de tests y declaradas en el propio código con su justificación: se relajaron después de un fallo, tras medir el ruido numérico que las explicaba (`tests/validate_lambda.py` y `residual_anatomy/diag_scaling_equivalence.py`).
3. **Detectar un resultado engañoso e identificar su origen**, en lugar de publicarlo o esconderlo.
4. **Trazabilidad completa**: cada resultado guardado con sus datos, su manifiesto y huellas criptográficas (SHA-256), para que un tercero pueda reproducirlo.
5. **Honestidad sobre el alcance**: separar lo PROBADO de la HIPÓTESIS y de lo ABIERTO, y conservar también los resultados negativos.

---

## 1. El problema, en términos sencillos

Un agujero negro atrapa la luz porque nada puede escapar más rápido que ella. Se puede crear algo parecido con sonido: si un gas fluye más rápido que la velocidad del sonido dentro de él, las ondas sonoras que intentan remontar la corriente quedan atrapadas. Hay un "horizonte" sónico.

En 1974, Stephen Hawking predijo que los agujeros negros emiten una débil radiación térmica. Su análogo sonoro debería emitir lo mismo: pares de vibraciones cuánticas, una escapa y la otra cae. En un gas de átomos enfriados casi al cero absoluto (un *condensado de Bose-Einstein*), ese efecto se puede calcular y, en principio, medir.

Este proyecto calcula esa radiación con un solver numérico propio. La magnitud clave es la **"temperatura"** de la radiación emitida, y la pregunta era cómo cambia cuando se modifica la física a escalas muy pequeñas. Es una versión de laboratorio de un problema abierto de la física teórica: el "problema transplanckiano".

## 2. Cómo se trabajó

**Tres papeles separados (más un asesor sin valor de revisión):**

- **Una persona (yo)** decide qué se hace, autoriza cada ejecución y arbitra los desacuerdos.
- **Claude Code (Anthropic)** implementa: escribe el código, lo ejecuta, verifica los números y redacta los informes.
- **ChatGPT (OpenAI)** actúa como **revisor externo**: revisa los diseños antes de que exista el código, audita el código antes de ejecutarlo y, al final, los paquetes de datos.
- **Claude en el chat de claude.ai (Anthropic)**, en la auditoría de noviembre, actuó como **asesor** de la persona: ayudó a coordinar y a decidir (por ejemplo, propuso guardar en un commit los criterios de E1-bis antes de ejecutarlo y bajar una tolerancia de un test). Es de la misma familia que el implementador, así que **sus sugerencias no cuentan como revisión independiente**; en este documento, «el revisor» es siempre ChatGPT.

El trabajo lo dirigió una persona sin formación académica en física. Por eso, el método se diseñó para que las conclusiones no dependieran de la autoridad de nadie, humana o artificial, sino de comprobaciones explícitas, documentadas y reproducibles.

**Reglas que se aplicaron sin excepciones:**

- **Diseño en papel → revisión → implementación → auditoría → ejecución.** Ningún cálculo importante se lanzó sin un diseño aprobado.
- **Criterios fijados de antemano** ("gates"). Antes de ejecutar se escribe qué resultado cuenta como éxito, cuál como fracaso y cuál como inconcluso. Si el resultado es desfavorable, se informa tal cual; no se mueve el umbral. (En los tests hubo dos excepciones declaradas: ver «Qué demuestra este proyecto», punto 2.)
- **Etiquetas obligatorias**: PROBADO / HIPÓTESIS / ABIERTO, y siempre "validado dentro del alcance documentado".
- **Persistencia**: cada resultado sale con sus CSV, un manifiesto con el entorno de software y huellas SHA-256 de código y datos.
- **Independencia real**: dos modelos de la misma familia no cuentan como dos verificaciones independientes, porque pueden compartir sesgos. Por eso el revisor es de otra empresa.

En total, el proyecto acumuló una veintena de gates documentados.

## 3. El resultado que parecía un descubrimiento

El simulador básico (versión 1.0) se validó primero y quedó congelado. Reproduce la predicción teórica de Finazzi y Parentani (2011): la temperatura obtenida es entre el 95,7 % y el 98,5 % de la temperatura de Hawking ideal, y se acerca al 100 % a medida que el horizonte es más suave. (Con la temperatura corregida que proponen esos autores, el acuerdo es mejor que el 0,75 %, pero esa cifra depende de un parámetro ajustado por ellos para sus propios perfiles.)

Después se añadió un parámetro, λ, con la intención de cambiar solo la física a escalas muy pequeñas sin tocar nada más. El resultado fue sorprendente: al alejar λ de su valor original, en cualquiera de los dos sentidos, **la temperatura de la radiación caía a cero**. El efecto esperado era de unos pocos por ciento; aquello era un "colapso" total.

Un resultado así podría haberse presentado como hallazgo. En cambio, se investigó (sección 4) y resultó ser un error del modelo. Y la investigación dejó una segunda lección, menos visible pero igual de importante: **una vez corregido el error, λ ya no hacía lo que se pretendía.** Se demostró que cambiar λ equivale exactamente a hacer el horizonte más ancho o más estrecho: los resultados con λ coinciden con los del simulador original en un horizonte reescalado, con una diferencia relativa de apenas unas milmillonésimas (3,5·10⁻⁹). Es decir, λ no aísla la física a escalas pequeñas; es otra forma de mover un parámetro que ya se conocía. La pregunta original, tal como estaba planteada, no tenía el contenido que parecía tener.

## 4. Cómo se identificó el error

1. **Pausa y auditoría bibliográfica.** Antes de construir un solver más complejo para estudiar el "colapso" (varios pasos de trabajo más), se decidió comprobar si la pregunta ya estaba respondida. Leyendo los artículos originales (Macher y Parentani 2009; Finazzi y Parentani 2012), resultó que **cálculos muy parecidos, con otros perfiles de flujo, ya estaban hechos y no mostraban ningún colapso**; nuestro resultado estaba en fuerte tensión con ellos.
2. **Hipótesis concreta.** Al cambiar λ, el código modificaba la ecuación de las vibraciones pero no actualizaba un término del estado de fondo del gas. El gas simulado ya no estaba en equilibrio.
3. **Test de simetría.** Si el fondo está en equilibrio, existe una vibración concreta, el *modo de Goldstone*, que debe tener frecuencia exactamente cero. Se midió:
   - el simulador original pasaba el test (el residuo bajaba al refinar la malla, como corresponde a un error de discretización);
   - la versión con λ lo violaba, y el residuo coincidía con la fórmula predicha en papel con un cociente de **1,0000**;
   - corrigiendo el término, el modo cero volvía.
4. **La prueba decisiva.** Se repitió el cálculo original con el fondo corregido, en las mismas condiciones. **El "colapso" desapareció.** La temperatura pasó a variar de forma suave, de 0,98 a 0,89 veces el valor de referencia, y coincide exactamente con el simulador original aplicado a un horizonte de anchura reescalada (sección 3).
5. **Consecuencia adicional.** El solver más complejo que estaba planificado (varios pasos de trabajo) **tenía el mismo defecto de diseño**. Se detectó antes de construirlo.

## 5. Validación externa

Que un instrumento coincida consigo mismo no basta: podría arrastrar un error sistemático. Se buscó un resultado **exacto** de otro grupo, para el mismo tipo de sistema: la configuración *waterfall* ("cascada") de Larré, Recati, Carusotto y Pavloff (2012): un flujo que cae por un escalón de potencial y cuya solución exacta se conoce, con fórmulas cerradas para la radiación.

Con criterios fijados antes de ejecutar, el simulador reprodujo:

- el coeficiente principal de la radiación con un error relativo de **3·10⁻⁵**;
- el coeficiente secundario con un error **inferior a 10⁻³**;

en los dos casos ensayados, con convergencia numérica del orden esperado y la sensibilidad al método de ajuste cuantificada.

Una auditoría posterior mostró que ese 3·10⁻⁵ venía casi entero de cómo se extrapolaba a frecuencia cero, no del simulador. Para comprobarlo sin ajustar nada a posteriori, se escribieron los criterios de una segunda prueba (E1-bis) y se guardaron en un commit **antes** de ejecutarla; ese commit es público (`001a7ac`). En dos casos que no se habían mirado, el simulador reprodujo el coeficiente principal con errores de **5·10⁻⁹ y 4·10⁻⁸**, y el secundario con **8·10⁻⁷ y 3·10⁻⁶**.

Las ecuaciones tomadas de los artículos se cotejaron después con su texto literal. El cotejo confirmó las convenciones del benchmark (la escala Λ de Finazzi y Parentani) y encontró una cita con el número de ecuación equivocado, ya corregida.

## 6. Lo que falló por el camino

Es la parte más reveladora del método. Siete errores reales, elegidos entre los veintiséis registrados (la lista completa está en el anexo):

| Qué falló | Quién lo detectó | Qué enseñó / qué cambió en el proceso |
|---|---|---|
| La temperatura teórica de referencia tenía un factor ½ de más; parecía una discrepancia física del 120 % | El revisor exigió auditar las convenciones; el implementador encontró el error contrastando con la literatura | Antes de interpretar una discrepancia como física, se auditan normalizaciones y convenciones |
| Una falsa alarma ("el condicionamiento explota") venía de un fallo de una librería numérica (OpenBLAS) al resolver el sistema adjunto | Auditoría adversaria en tres fases (implementador, segundo modelo, revisor externo) | Un número anómalo no basta para concluir: se triangula con una vía independiente (aquí, tres métodos certificados) |
| Un `max()` ocultaba valores NaN en el código de un gate (fue una de las razones de que ese gate necesitara diez versiones) | El revisor, auditando el código antes de ejecutar | Código "a prueba de fallos": cobertura obligatoria y NaN tratados de forma explícita |
| El "colapso" de la radiación (sección 3) | Auditoría bibliográfica pedida por la persona; el implementador identificó el defecto | Cualquier modelo modificado pasa un test de simetría antes de interpretar sus resultados |
| La primera medida de ese test estaba contaminada por los bordes del dominio | El implementador, al ver que el error crecía con la resolución cuando debía bajar | Se mide donde tiene sentido y la versión errónea se conserva como registro |
| Un control de convergencia daba cero por construcción: no comprobaba nada | El revisor, auditando el paquete de datos | Los controles también se auditan: órdenes calculados por separado y extrapolación declarada |
| Un control de canales no comprobaba lo que decía comprobar | El revisor, auditando el paquete | Comprobación explícita, registrada frecuencia a frecuencia |

Ninguno de estos errores llegó a una conclusión final. Los dos últimos afectaban a los controles, no a la física: tras corregirlos, los datos recalculados fueron idénticos byte a byte.

## 7. Fase 2: un intento de aislar la física de escalas pequeñas, sin veredicto

La sección 3 dejó una pregunta abierta: λ no aísla la física de escalas pequeñas porque solo mueve una escala. La fase 2 lo intentó con una dispersión de **dos escalas**, (ω − vk)² = c²k² + k⁴/4 + γc²k⁶, cuyo segundo número adimensional no lo absorbe ningún reescalado. Es un modelo fenomenológico de dispersión (en la línea de Corley y Jacobson), no una interacción física.

- **El instrumento quedó validado.** Una primera formulación falló sus gates por un problema de redondeo, y así consta. La versión revisada pasó los cinco: modo de Goldstone exacto, convergencia de orden 4, coincidencia con el simulador original a 6·10⁻⁹ cuando el término nuevo se anula, y no degeneración frente a cualquier reescalado.
- **La medición se pre-registró.** La predicción, las tolerancias y un **control** se pusieron en commit antes de calcular ninguna temperatura. La predicción era que el exponente del horizonte ensanchado pasara de 1/3 a 1/5, según un argumento de escalas que generaliza el de Coutant y Parentani. Tres addenda declararon, antes de medir, lo que se iba viendo, sin cambiar ningún criterio, y el revisor externo lo revisó todo antes de ejecutar.
- **Resultado: sin veredicto.** El control exigía que, sin el término nuevo, el método diera el exponente conocido (≈ 1/3). Dio 0,62. El observable elegido, la anchura de la ventana que reproduce la temperatura, mezcla la geometría del perfil con el ensanchamiento dispersivo, así que no medía lo que se pretendía. Como el control falló, la predicción no se evaluó.
- **Por qué importa:** con el término nuevo, el exponente salió 0,51, «en la dirección predicha». Sin un control fijado de antemano, habría sido tentador presentarlo como una confirmación. El control detuvo esa conclusión antes de que existiera. Es el mismo tipo de error que el de la fase 1, pero esta vez lo paró el propio método.
- *Datos descriptivos, sin conclusión:* con γ = 10, el término k⁶ baja la temperatura un 1,8 % en el horizonte más abrupto (w = 2) y un 0,04 % en el más suave (w = 8).

## 8. Qué no se afirma

- **No hay física nueva.** La robustez de la radiación frente a cambios a escalas pequeñas ya estaba establecida (Macher, Finazzi, Parentani y otros). El proyecto la reproduce; no la descubre.
- Que el "colapso" era un artefacto está demostrado **dentro del alcance ensayado** (un perfil de flujo, un rango de λ). Con el fondo corregido, los resultados coinciden con el simulador original reescalado hasta el redondeo, de modo que no queda margen para otro mecanismo en el alcance ensayado; la convergencia respecto al tamaño del dominio sigue sin estudiarse.
- El parámetro λ no estudia la física a escalas pequeñas con la geometría fija: con el fondo corregido es un reescalado de la anchura del horizonte. Las conclusiones sobre correlaciones obtenidas con λ se han retirado.
- La validación externa cubre el simulador básico (λ = 1), en la configuración y los dos casos ensayados.
- Algunos resultados intermedios quedaron como inconclusos o negativos, y así constan. Por ejemplo, un gate de referencia de fase terminó en "clase C" (criterio no superado) y se archivó tal cual.
- La fase 2 no da ningún resultado sobre la física de escalas pequeñas: terminó sin veredicto.

## 9. En una frase

Una persona organizó a dos IAs con papeles separados y reglas fijadas de antemano. El sistema produjo un resultado engañoso, lo detectó, identificó su origen, lo corrigió y validó el instrumento contra un resultado exacto publicado, con cada paso reproducible. En una segunda fase, un control fijado de antemano impidió sacar una conclusión que los datos no justificaban.

---

*Herramientas: implementación con Claude Code (Anthropic); revisión externa con ChatGPT (OpenAI); asesoría de coordinación con Claude en claude.ai (Anthropic), sin valor de revisión independiente. Código propio en Python (NumPy, SciPy).*

---

## Anexo técnico (para físicos)

**Sistema.** Condensado de Bose-Einstein 1D, ecuaciones de Bogoliubov-de Gennes (BdG) estacionarias; ℏ = m = g = 1. Perfil de densidad n(x) en tanh (Mach 0,5 aguas arriba, 2 aguas abajo, anchura w = 2); no confundir con la configuración *waterfall* de Larré et al. usada en E1. Matriz S 3×3 por *matching* modal exacto con operador en diferencias finitas de segundo orden. Observable: n_L(ω), la ocupación del canal Hawking aguas arriba (Finazzi y Parentani 2011, Ec. 42). T_eff = lím ω·n_L.

**Benchmark v1.0 (congelado).** T_up/(κ_local/2π) = 0,957 → 0,985 para w ∈ {2, 3, 4, 6, 8}, frente a Finazzi y Parentani, PRD 83, 084010 (2011). Con su κ̄ de horizonte ensanchado, T_up/(κ̄/2π) ∈ [0,9926, 1,0050]; esta segunda cifra es condicional a su anchura de ventana d = 1,58 (con d = 1 o 2, la desviación máxima es del 2,4 % y del 3,4 %). La escala Λ = 2mc²/ℏ de su Ec. 9 se ha cotejado con el artículo.

**El artefacto.** El operador con dispersión (ω−vk)² = c²k² + λk⁴/4 equivale a BdG con masa m′ = 1/√λ y acoplo g′n = c²/√λ, pero mantenía fija la presión cuántica del fondo, Q_p, en lugar de Q_p,λ = √λ·Q_p. Residuo de Goldstone sobre Φ_G = (√n, −√n): ½(1 − √λ)(√n)″. Medido / predicho = 1,0000, independiente de la malla. Con Q_p,λ el residuo baja a O(dx²). El operador Padé planificado tenía el residuo (K_P − K_1)√n (cociente medido ≈ 1,00).

**T2: n_L con el fondo corregido** (w = 2, ω/(κ/2π) ∈ [0,05, 0,27], N = 9600 y 14400):

| λ | 0,25 | 0,5 | 1 | 2 | 4 |
|---|---|---|---|---|---|
| T_eff/T_H, fondo inconsistente | −0,006 | −0,027 | 0,957 | −0,053 | −0,017 |
| T_eff/T_H, fondo consistente | 0,979 | 0,971 | 0,957 | 0,932 | 0,892 |

(T_H = κ_local/2π. Con el fondo inconsistente, n_L ∝ ω a baja frecuencia, de modo que el T_eff extrapolado sale ≈ 0, ligeramente negativo por la extrapolación. Con el fondo consistente, la diferencia entre mallas es ≤ 6·10⁻⁵.)

**Equivalencia λ ↔ anchura.** Con el fondo consistente, x′ = x/√λ, ω′ = √λ·ω transforman el operador con λ sobre un perfil de anchura w en el operador físico (λ = 1) sobre un perfil de anchura w/√λ: M₁(x′) = √λ·M_λ(x), también en la malla (L′ = L/√λ, mismo N). Comprobado en los 90 puntos de T2: diferencia relativa en n_L ≤ 3,2·10⁻⁶ con el solve estándar (su ruido de redondeo) y ≤ 3,5·10⁻⁹ con equilibrado y refinamiento iterativo. La fila «fondo consistente» es, por tanto, la curva T/T_H frente a la anchura w/√λ = 4 → 1.

**E1: validación frente a Larré et al., PRA 85, 013621 (2012).** a = |f_{u,d2}|² (su Ec. A8) y b = 2Re(f*h) (Ec. A14), extraídos de ω·n_L = a + bω + cω² + dω³:

| m_u | a numérico / exacto | error | b numérico / exacto | error |
|---|---|---|---|---|
| 0,5 | 0,131737 / 0,131734 | 2,4·10⁻⁵ | −0,490187 / −0,489796 | 8,0·10⁻⁴ |
| 0,7 | 0,066910 / 0,066908 | 2,9·10⁻⁵ | −0,498434 / −0,497974 | 9,2·10⁻⁴ |

Γ y T_H (su Ec. 57) reproducidos. Integridad: regresión del driver frente a v1.0 = 0; residuo Gross-Pitaevskii de orden 2; defecto de pseudo-unitariedad con orden observado entre 1,87 y 1,98; canales verificados en 240/240 frecuencias.

**E1-bis (criterios pre-registrados en el commit `001a7ac`, anterior a la ejecución).** El error de E1 resultó ser el sesgo de la extrapolación ω → 0 (ventana W1 con ajuste cúbico), no redondeo (10⁻⁹–10⁻⁷ en esa malla) ni discretización. Ventana W2 = [0,05; 0,3]·T_H con ajuste de orden 4, dx = 0,01/0,005/0,0025, solve refinado y Richardson con p = 2:

| m_u | error en a (A8) | error en b (A14) | orden observado |
|---|---|---|---|
| 0,6 | 5,4·10⁻⁹ | 7,9·10⁻⁷ | 2,04 |
| 0,7 | 4,2·10⁻⁸ | 3,1·10⁻⁶ | 2,00 |

Criterio principal: error en a ≤ 10⁻⁶ en las dos, y control de ventana (W2 frente a W3) ≤ 10⁻⁶: PASS. m_u = 0,5 se repitió como réplica, sin contar para el veredicto (1,1·10⁻⁸).

**Nota sobre las dos temperaturas de referencia.** En el benchmark v1.0 el cociente es frente a T_κ̄ = κ̄/(2π), la temperatura con el gradiente promediado de Finazzi y Parentani (2011); por eso sale ≈ 1,00. En T2 el cociente es frente a T_local = κ_local/(2π), con el gradiente en el horizonte; por eso λ = 1 da 0,957. No son resultados contradictorios: son cocientes con denominadores diferentes.

**Registro completo de errores (14).** Los cuatro últimos afectaban a los controles y a la documentación, no a la física: tras corregirlos, los datos recalculados fueron idénticos byte a byte.

| Qué falló | Quién lo detectó | Qué cambió en el proceso |
|---|---|---|
| La temperatura teórica de referencia tenía un factor ½ de más; parecía una discrepancia física del 120 % | El revisor exigió auditar las convenciones; el implementador encontró el error contrastando con la literatura | Regla: auditar siempre normalizaciones y convenciones antes de interpretar una discrepancia |
| Un ajuste con ventana de frecuencias fija falló al cambiar la escala del problema | Los residuos del ajuste (salvaguarda pedida por el revisor) | Ventanas de ajuste escaladas con la física del caso |
| Una contabilidad de flujo usaba el signo equivocado (déficit del 2 %) | El revisor indicó el paso a comprobar; el implementador lo corrigió | Distinguir siempre magnitudes que parecen intercambiables y no lo son |
| Dos falsas alarmas: una "explosión" del condicionamiento y una no-convergencia | Auditoría adversaria en tres fases (implementador, segundo modelo, revisor externo) | Triangular con una vía independiente antes de concluir nada de un número anómalo |
| El defecto de fondo en una de ellas: un fallo de una librería numérica (OpenBLAS) al resolver el sistema adjunto | Esa misma auditoría | Ruta numérica alternativa certificada con tres métodos independientes |
| El solver "estabilizado" validado era de una versión simplificada del modelo, no de la que hacía falta | Auditoría adversaria | Validar el instrumento que de verdad se va a usar, no el más fácil |
| Resultados que solo se imprimían en pantalla y se perdían | Auditoría adversaria | Todo resultado se guarda con CSV, manifiesto y SHA-256 |
| Fallos de robustez en el código de un gate (diez versiones antes de aprobarlo): por ejemplo, un `max()` que ocultaba valores NaN | El revisor, auditando el código antes de ejecutar | Código "a prueba de fallos": cobertura obligatoria y NaN tratados explícitamente |
| El "colapso" de la radiación (sección 3) | Auditoría bibliográfica pedida por la persona; el implementador identificó el defecto | Test de simetría obligatorio para cualquier modelo modificado |
| La primera medida del test de simetría estaba contaminada por los bordes del dominio | El implementador, al ver que el error crecía con la resolución cuando debía bajar | Medir en el interior y conservar la versión errónea como registro |
| Un control de convergencia daba cero por construcción: no comprobaba nada | El revisor, auditando el paquete de datos | Órdenes de convergencia calculados por separado y extrapolación declarada |
| Un control de canales no comprobaba lo que decía comprobar | El revisor, auditando el paquete | Comprobación explícita, registrada frecuencia a frecuencia |
| El paquete entregado no se podía reproducir solo | El revisor | Paquetes autocontenidos, con dependencias y sumas de verificación |
| Una huella criptográfica atribuida al nombre de archivo equivocado | El revisor | Separar la versión ejecutada de la versión actualizada de cada documento |

**Errores encontrados en la auditoría de noviembre de 2026 (5).** El segundo retiró conclusiones de documentos del repositorio (no de este texto); los demás no cambian ningún número publicado.

| Qué falló | Quién lo detectó | Qué cambió en el proceso |
|---|---|---|
| Los tests imprimían resultados pero nunca fallaban: no tenían aserciones | El implementador, en la auditoría de solo lectura pedida por un asesor | Todos los tests comprueban con tolerancia y devuelven error; `tests/run_all.py` los ejecuta todos |
| Las correlaciones calculadas con λ heredaban el mismo defecto de fondo que el "colapso" | El implementador, en esa misma auditoría, al repasar qué resultados usaban el operador defectuoso | Esas conclusiones se retiraron, con aviso público |
| La amplitud absoluta de la correlación densidad-densidad salía a la mitad: el prefactor del artículo corresponde a otra normalización de modos. Además, la ecuación se citaba con un número equivocado | El implementador marcó el posible factor 4π en la auditoría de la fase 1; el revisor localizó la fórmula (Ec. 25) y el número de ecuación equivocado; el implementador lo zanjó con una prueba exacta en un condensado homogéneo | Test de normalización absoluta con tolerancia; ecuaciones cotejadas con el texto literal de los artículos |
| Una figura publicada se había generado con una versión del código anterior a la que describía el texto (sin la ventana de suavizado) | El implementador, al regenerarla | Las figuras se regeneran con los scripts del repositorio antes de publicar |
| Una nota de auditoría atribuía al cálculo E1 un ruido de redondeo medido en otra configuración | El implementador, al medirlo en la configuración correcta | Medir cada magnitud en la configuración a la que se refiere |

**Errores de la fase 2 (7).** Ninguno llegó a una conclusión publicada.

| Qué falló | Quién lo detectó | Qué cambió en el proceso |
|---|---|---|
| El *matching* del primer diseño dejaba el sistema casi singular | El implementador, en las pruebas con flujo uniforme, antes de los gates | Contar el estado del sistema discreto antes de fijar las condiciones de contorno |
| La primera formulación tenía un suelo de redondeo que crecía al refinar la malla (∝ dx⁻⁴) | Los gates F0 y F1 | Reformulación (variables s, d; cuarto orden; mallas más gruesas) con criterios nuevos en commit antes de implementarla |
| Los coeficientes de las diferencias finitas se creaban con 15 dígitos en vez de 60 | El gate F1′, en su primera ejecución | Registrar el fallo, corregir el error y repetir con los mismos criterios |
| El modelo se presentó como una interacción física | El asesor | Es un modelo fenomenológico de dispersión; corregido en el diseño |
| Se dio por «dominado por k⁶» un caso que estaba en la transición | Pregunta del asesor; cálculo del implementador | Comprobar el régimen en papel antes de medir (addendum 2) |
| Faltaba un término en la expansión de la dispersión | El revisor externo | Análisis de sensibilidad (addendum 3); no cambiaba los casos con k⁶ |
| El observable no aislaba el ensanchamiento dispersivo | El control K pre-registrado | La medición se declaró sin veredicto |

**Referencias.**
- J. Macher, R. Parentani, PRA 80, 043601 (2009)
- A. Recati, N. Pavloff, I. Carusotto, PRA 80, 043603 (2009)
- S. Finazzi, R. Parentani, PRD 83, 084010 (2011)
- S. Finazzi, R. Parentani, PRD 85, 124027 (2012)
- A. Coutant, R. Parentani, S. Finazzi, PRD 85, 024021 (2012)
- P.-É. Larré, A. Recati, I. Carusotto, N. Pavloff, PRA 85, 013621 (2012)

**Dónde está cada cosa** (en el repositorio del proyecto, <https://github.com/jorditorresribas/bec-analog-hawking>):
- `docs/FP_BENCHMARK.md`: benchmark v1.0.
- `residual_anatomy/NOVELTY_AUDIT.md`: auditoría de novedad, con el Apéndice A (derivación, test de Goldstone, T2).
- `residual_anatomy/E1_RESULT.md` y `residual_anatomy/GATE_EXT_LRCP_DESIGN.md`: validación externa; `residual_anatomy/E1BIS_PREREG.md`: criterios de E1-bis.
- `docs/CORRELATIONS.md` §9: equivalencia λ ↔ anchura.
- `tests/run_all.py`: todas las validaciones.
- `docs/FASE2_DISENO.md`, `residual_anatomy/FASE2_PREREG.md` (y addenda), `residual_anatomy/FASE2_GATES_RESULT.md` y `residual_anatomy/FASE2_RESULT.md`: fase 2.
- `residual_anatomy/P2_RESULT.md`: un resultado negativo archivado.
- `out/`: CSV y manifiestos.
